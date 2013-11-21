/*
 * Copyright 2012-2013 TopCoder, Inc.

 *
 * This code was developed under U.S. government contract NNH10CD71C. 

 *
 * Licensed under the Apache License, Version 2.0 (the "License");

 * You may not use this file except in compliance with the License.

 * You may obtain a copy of the License at:
 *     http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software

 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

 * See the License for the specific language governing permissions and

 * limitations under the License.
 */

#include <vector>
#include <iostream>
#include <sstream>
#include <string> 
#include <algorithm>
#include <map>
#include <memory.h>

using namespace std;

namespace place2 {

//#define DEBUG
#ifdef DEBUG
#define dout(x) cout<<x
#endif
#ifndef DEBUG
#define dout(x)
#endif


double PARAMA = 0.5; // Potenz
double PARAMB = 2.5; // Multiplikator
double PARAMC = 60; // Smooth day
double PARAMD = 180.0; // Param exp abhahme Teile 
double PARAME = 30.0; // Param exp abnahme Minus days
double PARAMF = 0.25; // Faktor social rate
double PARAMG = 0.65; // Influence social rate
double countryInf = 0.03;
double neighborsInf = 0.16;


/*
double PARAMA; // Potenz
double PARAMB; // Multiplikator
double PARAMC; // Smooth day
double PARAMD; // Param exp abhahme Teile 
double PARAME; // Param exp abnahme Minus days
double PARAMF; // Faktor social rate
double PARAMG; // Influence social rate
double countryInf = 0.03;
double neighborsInf = 0.16;
*/


int toInt(string input) {
	if( input == "_" ) {
		return -1;
	}
	stringstream stream(input);
	int result;
	stream>>result;
	return result;
}


double toDouble(string input) {
	stringstream stream(input);
	double result;
	stream>>result;
	return result;
}

void trim(string& str) {
	str.erase(str.find_last_not_of(" \n\r\t")+1);
}

void split(vector<string> &tokens, const string &text, char sep) {
	tokens.clear();
	size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != string::npos) {
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
}



const int COUNTRIES = 254;
const int REGIONS = 3671;
const int MAX_DAYS = 10000;


class Coord {
public:
	double x,y;
	Coord() :x(0), y(0) {}
	Coord(double x_, double y_): x(x_), y(y_) {}
	
	Coord(string str) {
		if( str == "_" ) {
			x = NAN;   // nan 
			y = NAN;   // nan 
			return;
		}
		vector<string> tokens;
		split(tokens, str, ',');
		x = toDouble(tokens[0]);
		y = toDouble(tokens[1]);
	}
	
	Coord operator+(Coord obj) {
		return Coord(x + obj.x, y + obj.y);
	}
	Coord operator*(double skalar) {
		return Coord(x*skalar, y*skalar);
	}
	
	double operator|(Coord obj) {
		return (x-obj.x)*(x-obj.x) + (y-obj.y)*(y-obj.y);
	}
};

class Polygon {
	public:
		vector<Coord> mCoords;
		Coord mCenter;
		double mArea;
		void calcArea() {
			mArea = 0.0;
			for( size_t i = 0; i < mCoords.size()-1; i++) {
				mArea += ( mCoords[i+1].x*mCoords[i].y - mCoords[i].x * mCoords[i+1].y);
			}
			mArea *= 0.5;
		}
		void calcCenter() {
			mCenter = Coord(0,0);
			for( size_t i = 0; i < mCoords.size()-1; i++) {
				mCenter = mCenter + (mCoords[i] + mCoords[i+1])*(mCoords[i+1].x * mCoords[i].y - mCoords[i].x*mCoords[i+1].y);
			}
			mCenter =  mCenter * ( 1.0 / (6.0 * mArea));
		}
		Polygon(vector<Coord>& coords) {
			mCoords = coords;
			calcArea();
			calcCenter();
		}
		
		const double EPS =  0.00001;
		bool shareVertex(Polygon& obj) {
			double dis = 0.0;
			for( Coord& a: mCoords) {
				for( Coord& b: obj.mCoords) {
					dis = sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
					if( dis < EPS ) {
						return true;
					}
				}
			}
			return false;
		}
		
		Coord getCenter() {	return mCenter; }
		double getArea() {	return mArea;	}
};

int dayID;
int firstDayID;

class Location {
public:
	int mPrecision;
	Coord mCoord;
	int mCountryId;
	int mRegionId;
	Location() {};
	Location( int precision, Coord coord, int countryId, int regionId):
				mPrecision(precision), mCoord(coord), mCountryId(countryId), mRegionId(regionId) {}
};

class Participant{
public:
	int mId;
	Location mLocation;
	Participant() {};
	Participant(int id, Location location): mId(id), mLocation(location) {}
};map<int,Participant> participants;

class SocialEvent {
public:
	int mP1,mP2;
	int mType;
	Location mLocation;
	bool mImportant;
	int mMediaCoverage;
	int mMediaSentiments;
	SocialEvent(int p1, int p2, int type, Location location, bool important, int mediaC, int mediaS): 
				mP1(max(p1,p2)), mP2(min(p1,p2)), mType(type), mLocation(location), mImportant(important), mMediaCoverage(mediaC), mMediaSentiments(mediaS) {}
};vector<SocialEvent> socialEvents;

class Segment {
	public:
		vector<int> mEventDays;
		int mSocialEventCount[MAX_DAYS];
		int mSocialEventTotal;
		
		Segment() {
			for( int i = 0; i < MAX_DAYS; i++) {
				mSocialEventCount[i] = 0;
			}
			mSocialEventTotal = 0;
		}
		
		virtual void registerEvent() {
			mEventDays.push_back(dayID);
		}
		
		virtual void registerSocialEvent() {
			mSocialEventCount[dayID - firstDayID]++;
			mSocialEventTotal++;
		}
		
		double getWGH() {
			if( mEventDays.size() == 0 ) {
				return 0;
			}
			return 1.0 - tanh((dayID - mEventDays[mEventDays.size()-1] + 10.0) / 180.0);
		}
		
		double getWGH_() {
			if( mEventDays.size() == 0 ) {
				return 1.0;
			}
			return tanh((dayID - mEventDays[mEventDays.size()-1] + 10.0) / 180.0);
		}
		
		double getExpAbnahme() {
			if( mEventDays.size() == 0 ) {
				return 0.0;
			}
			double d = dayID - mEventDays[mEventDays.size()-1] - PARAME;
			return pow(2.0, -d/PARAMD) * 1.35;
		}
		
		
		double getMaxEventProb() {
			double res = 0.0;
			for( size_t i = 0; i < mEventDays.size(); i++) {
				double nenner =  PARAMC + dayID - mEventDays[i];
				res = max(res, (double)(mEventDays.size() - i - 1) / nenner);
			}
			return pow(res, PARAMA) * PARAMB;
		}
	
		
		static const int MEDIA_TIME = 5;
		double getSocialEventRate() {
			double sumTime = 0;
			for( int  i = 0; i < MEDIA_TIME; i++) {
				sumTime += (double)mSocialEventCount[dayID - firstDayID - i];
			}
			if( mSocialEventTotal == 0 || sumTime == 0 ) {
				return 0.0;
			}
			sumTime /= (double)MEDIA_TIME;
			double average = (double)mSocialEventTotal / (double)(dayID - firstDayID);
			//if( sumTime > 0.01)
			//	cout<<sumTime/average - 3.0<<endl;
			return sumTime / average;
		}
		
		double getPrediction() {
			double rate = getSocialEventRate();
			double eventProb = getMaxEventProb();
			//rate /= 8.0;
			eventProb *= ((1.0 - PARAMG) + PARAMG*tanh(rate*PARAMF) );
			
			/*
			if( rate > 12.0 ) {
				eventProb =  pow(eventProb, 0.625) ;
			}else if( rate > 9.0 ) {
				eventProb = pow(eventProb, 0.725) ;
			}else if( rate > 6.0 ) {
				eventProb = pow(eventProb, 0.825) ;
			}else if( rate > 3.0 ) {
				eventProb = pow(eventProb, 0.925);
			}*/
			eventProb *= getExpAbnahme();
			return min(1.0, max(0.0, eventProb ));
		}
			
};


class Country : public Segment{
	public: 
		int mId;
		string mName;
		vector<int> mRegions;
		
		Country(int id, string name): Segment(), mId(id), mName(name) {}
		
};vector<Country> countries;

class Region : public Segment{
	private:
		int mCountryId;
		string mCountryName;
		int mID;
		string mName;
		
		
		
		vector<Polygon> mInner;
	public:
		vector<Polygon> mOuter;
		vector<int> mNeighbors;
		Coord mCenter;
		double mArea;

		Region(	int countryID, string countryName, int id, string name, vector<Polygon>& inner, vector<Polygon>& outer):
				Segment(), mCountryId(countryID), mCountryName(countryName), mID(id), mInner(inner), mOuter(outer) {
			mCenter = Coord(0,0);
			mArea = 0.0;
			for( Polygon& p: outer ) {
				mArea += p.mArea;
				mCenter = mCenter + p.mCenter * mArea;
			}
			if( mArea < 0.000001 ) {
				mArea = 0.00001;
				cout<<"Bad\n";
			}
			mCenter = mCenter * ( 1.0 / mArea);
		}
		
		void print() {
			dout(mCountryName<<" "<<mName<<endl);
			dout("Center:"<<mCenter.x<<","<<mCenter.y<<endl);
			dout("Area:"<<mArea<<endl);
		}
		
		
		void registerEvent() {
			Segment::registerEvent();
			countries[mCountryId].registerEvent();
		}
		
		void registerSocialEvent() {
			Segment::registerSocialEvent();
			countries[mCountryId].registerSocialEvent();
		}
		
		int getCountryID() {	return mCountryId; }
		
};vector<Region> regions;


double MAX_DIS = 1.2;
class MassAtrocityPredictor : public BaseMassAtrocityPredictor {
private:
	void calcNeighbords() {
		for( int i = 0; i < REGIONS; i++) {
			for( int j = i+1; j < REGIONS; j++) {
				/*
				Coord vec = regions[i].mCenter + (regions[j].mCenter*(-1.0));
				double dis = sqrt(vec.x*vec.x +  vec.y*vec.y);
				if( dis < MAX_DIS ) {
					regions[i].mNeighbors.push_back(j);
					regions[j].mNeighbors.push_back(i);
				}
				*/
				Coord vec = regions[i].mCenter + (regions[j].mCenter*(-1.0));
				double dis = sqrt(vec.x*vec.x +  vec.y*vec.y);
					if( dis < MAX_DIS*5 ) {
					for( Polygon& p1: regions[i].mOuter ) {
						for( Polygon& p2: regions[j].mOuter ) {
							if( p1.shareVertex(p2) ) {
								regions[i].mNeighbors.push_back(j);
								regions[j].mNeighbors.push_back(i);
							}
						}
					}
				}
			}
		}
	}
	
	void loadRegions(vector<string>& data) {
		vector<string> tokens;
		for( size_t i = 0; i < data.size();) {
			if( data[i] != "" ) {
				// Country 
				split(tokens, data[i], ' ');
				int cid = toInt(tokens[0]);
				string cname = tokens[1];
				if( cid >= (int)countries.size() ) {
					countries.push_back(Country(cid, cname));
				}
				i++;
				// Region 
				split(tokens, data[i], ' ');
				int rid = toInt(tokens[0]);
				string rname = tokens[1];
				i++;
				// Polygons
				vector<Polygon> inner,outer;
				while( i < data.size() ) {
					bool binner;
					string instr = data[i];
					trim(instr);
					if( instr == "inner" ) {
						binner = true;
					}else if( instr == "outer" ) {
						binner = false;
					}else{
						break;
					}
					i++;
					split(tokens, data[i], ' ');
					vector<Coord> coords(tokens.size() + 1);
					vector<string> tmp;
					for( size_t i = 0; i < tokens.size(); i++) {
						split(tmp, tokens[i], ',');
						coords[i] = Coord(toDouble(tmp[0]), toDouble(tmp[1]));
					}
					coords[coords.size()-1] = coords[0];
					if( binner ) {
						inner.push_back(Polygon(coords));
					}else {
						outer.push_back(Polygon(coords));
					}
					i++;
				}
				regions.push_back(Region(cid, cname, rid, rname, inner, outer));
			}else {
				i++;
			}
		}
		calcNeighbords();
	}
	
	
	
	void handleEvents(vector<string>& data) {
		for( string line: data) {
			vector<string> tokens;
			split(tokens, line, ' ');
			int regionID = toInt(tokens[3]);
			regions[regionID].registerEvent();
		}
	}
	
	Location makeLocation(string& precision, string& locdesc, string& cid, string& rid) {
		return Location( toInt(precision), Coord(locdesc),  toInt(cid), toInt(rid));
	}
	
	void registerParticipant(int id, Location location) {
		if( participants.count(id) == 0 ) {
			participants[id] = Participant(id, location);
		}
	}
	
	void handleSocialEvents(vector<string>& data) {
		vector<string> tokens;
		for( string line: data) {
			split(tokens, line, ' ');
			
			switch( tokens[10][0] - 'a' ) {
				case 3:
				//case 4:
				//case 5:
				case 7:
				case 8:
				case 11:
				//case 12:
				//case 13:
				//case 18:
				//case 19:
				
					if( /* tokens[15][0] == 't' && */tokens[14][0] != '_' ) {
						/*
						// Participants
						int p1 = -1;
						int p2 = -1;
						if( tokens[0] != "_" ) {
							p1 = toInt(tokens[0]);
							registerParticipant(p1, makeLocation(tokens[1], tokens[2], tokens[3], tokens[4]));
						}
						if( tokens[5] != "_" ) {
							p2 = toInt(tokens[5]);
							registerParticipant(p2, makeLocation(tokens[6], tokens[7], tokens[8], tokens[9]));
						}
						int type = tokens[10][0]-'a';
						Location location = makeLocation(tokens[11], tokens[12], tokens[13], tokens[14]);
						bool important = (tokens[15][0] == 't');
						int mediaC = toInt(tokens[16]);
						int mediaS = toInt(tokens[17]);
						SocialEvent event = SocialEvent(p1, p2, type, location, important, mediaC, mediaS);
						
						regions[event.mLocation.mRegionId].registerSocialEvent(socialEvents.size());
						socialEvents.push_back(event);
						*/
						regions[toInt(tokens[14])].registerSocialEvent();
					}else if(tokens[13][0] != '_') {
						countries[toInt(tokens[13])].registerSocialEvent();
					}
				default:
					break;
			}
		}
	}
	
	double getPrediction(int regionID) {
		double regionP = regions[regionID].getPrediction();
		double countryP = countries[regions[regionID].getCountryID()].getPrediction();
		
		// Neighbors
		if( regions[regionID].mNeighbors.size() > 0 ) {
			double neighborsP = 0.0;
			for( int i: regions[regionID].mNeighbors ) {
				neighborsP += regions[i].getPrediction();
			}
			neighborsP /= (double)regions[regionID].mNeighbors.size();
			regionP = regionP * (1.0 - neighborsInf) + neighborsP*neighborsInf;
		}
		return (regionP * (1.0 - countryInf) + countryP * countryInf );
	}
	
public:
	int receiveData(int dataSourceId, int dayID_, vector<string> data) {
		dayID = dayID_;
		switch( dataSourceId ) {
			case 0: // Atrocities
				dout("HandleEvents..\n");
				handleEvents(data);
				dout(".done\n");
				break;
			case 1: // Sociopolitical
				dout("Handle Social events..\n");
				handleSocialEvents(data);
				dout(".done\n");
				break;
				
			case 2:// Regions
				dout("Load regions..\n");
				firstDayID = dayID;
				loadRegions(data);
				dout("..done\n");
				break;
			default:
				break;
		}
		return 0;
	}
	
	vector<double> predictAtrocities(int dayID_) {
		dout("Predict events..\n");
		dayID = dayID_;
		vector<double> result;
		for( int i = 0; i < REGIONS; i++) {
			result.push_back(getPrediction(i));
		}
		dout("done\n");
		return result;
	}
	
	void finish() {
		for( int d = firstDayID; d < dayID; d++) {
			for( Region r: regions ) {
				cout<<r.mSocialEventCount[d-firstDayID]<<";";
			}
			cout<<endl;
		}
	}
};

}
