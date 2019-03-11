% Calculate kinetic energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:10
% EndTime: 2019-03-08 21:40:14
% DurationCPUTime: 3.93s
% Computational Cost: add. (2507->345), mult. (5854->479), div. (0->0), fcn. (6982->10), ass. (0->158)
t371 = Icges(4,1) + Icges(5,2);
t370 = Icges(5,1) + Icges(4,3);
t369 = Icges(6,1) + Icges(7,1);
t368 = -Icges(4,4) - Icges(5,6);
t367 = Icges(5,4) - Icges(4,5);
t366 = Icges(6,4) + Icges(7,4);
t365 = Icges(5,5) - Icges(4,6);
t364 = Icges(6,5) + Icges(7,5);
t363 = Icges(4,2) + Icges(5,3);
t362 = Icges(6,2) + Icges(7,2);
t361 = Icges(6,6) + Icges(7,6);
t360 = Icges(6,3) + Icges(7,3);
t359 = rSges(7,3) + qJ(6);
t287 = sin(pkin(10));
t289 = cos(pkin(10));
t295 = cos(qJ(2));
t290 = cos(pkin(6));
t293 = sin(qJ(2));
t320 = t290 * t293;
t252 = t287 * t295 + t289 * t320;
t288 = sin(pkin(6));
t332 = cos(qJ(3));
t309 = t288 * t332;
t331 = sin(qJ(3));
t234 = t252 * t331 + t289 * t309;
t319 = t290 * t295;
t251 = t287 * t293 - t289 * t319;
t292 = sin(qJ(5));
t294 = cos(qJ(5));
t201 = t234 * t294 - t251 * t292;
t326 = t234 * t292;
t202 = t251 * t294 + t326;
t308 = t288 * t331;
t235 = t252 * t332 - t289 * t308;
t357 = t201 * t361 + t202 * t364 + t235 * t360;
t254 = -t287 * t320 + t289 * t295;
t236 = t254 * t331 - t287 * t309;
t253 = t287 * t319 + t289 * t293;
t203 = t236 * t294 - t253 * t292;
t325 = t236 * t292;
t204 = t253 * t294 + t325;
t237 = t254 * t332 + t287 * t308;
t356 = t203 * t361 + t204 * t364 + t237 * t360;
t355 = t201 * t362 + t202 * t366 + t235 * t361;
t354 = t203 * t362 + t204 * t366 + t237 * t361;
t353 = t201 * t366 + t202 * t369 + t235 * t364;
t352 = t203 * t366 + t204 * t369 + t237 * t364;
t351 = t234 * t363 + t235 * t368 + t251 * t365;
t350 = t236 * t363 + t237 * t368 + t253 * t365;
t349 = t234 * t365 - t235 * t367 + t251 * t370;
t348 = t236 * t365 - t237 * t367 + t253 * t370;
t347 = t368 * t234 + t235 * t371 - t367 * t251;
t346 = t368 * t236 + t237 * t371 - t367 * t253;
t258 = -t290 * t332 + t293 * t308;
t321 = t288 * t295;
t238 = t258 * t294 + t292 * t321;
t324 = t258 * t292;
t239 = -t294 * t321 + t324;
t259 = t290 * t331 + t293 * t309;
t345 = t238 * t361 + t239 * t364 + t259 * t360;
t344 = t238 * t362 + t239 * t366 + t259 * t361;
t343 = t238 * t366 + t239 * t369 + t259 * t364;
t342 = t258 * t363 + t259 * t368 - t321 * t365;
t341 = t368 * t258 + t259 * t371 + t367 * t321;
t340 = t258 * t365 - t259 * t367 - t321 * t370;
t330 = pkin(7) * t290;
t329 = pkin(5) * t294;
t327 = Icges(2,4) * t287;
t323 = t287 * t288;
t322 = t288 * t289;
t318 = rSges(7,1) * t202 + rSges(7,2) * t201 + pkin(5) * t326 + t235 * t359 + t251 * t329;
t317 = rSges(7,1) * t204 + rSges(7,2) * t203 + pkin(5) * t325 + t237 * t359 + t253 * t329;
t316 = rSges(7,1) * t239 + rSges(7,2) * t238 + pkin(5) * t324 + t259 * t359 - t321 * t329;
t315 = qJD(2) * t288;
t314 = V_base(5) * qJ(1) + V_base(1);
t310 = qJD(1) + V_base(3);
t267 = t287 * t315 + V_base(4);
t279 = qJD(2) * t290 + V_base(6);
t233 = qJD(3) * t253 + t267;
t266 = -t289 * t315 + V_base(5);
t232 = qJD(3) * t251 + t266;
t255 = -qJD(3) * t321 + t279;
t261 = pkin(1) * t287 - pkin(7) * t322;
t307 = -t261 * V_base(6) + V_base(5) * t330 + t314;
t262 = pkin(1) * t289 + pkin(7) * t323;
t306 = V_base(4) * t261 - t262 * V_base(5) + t310;
t305 = V_base(6) * t262 + V_base(2) + (-qJ(1) - t330) * V_base(4);
t225 = pkin(2) * t252 + pkin(8) * t251;
t260 = (pkin(2) * t293 - pkin(8) * t295) * t288;
t304 = -t225 * t279 + t266 * t260 + t307;
t226 = pkin(2) * t254 + pkin(8) * t253;
t303 = t267 * t225 - t226 * t266 + t306;
t302 = t279 * t226 - t260 * t267 + t305;
t227 = pkin(3) * t259 + qJ(4) * t258;
t301 = qJD(4) * t236 + t232 * t227 + t304;
t196 = pkin(3) * t235 + qJ(4) * t234;
t300 = qJD(4) * t258 + t233 * t196 + t303;
t197 = pkin(3) * t237 + qJ(4) * t236;
t299 = qJD(4) * t234 + t255 * t197 + t302;
t205 = pkin(4) * t251 + pkin(9) * t235;
t241 = -pkin(4) * t321 + t259 * pkin(9);
t298 = t232 * t241 + (-t196 - t205) * t255 + t301;
t206 = pkin(4) * t253 + pkin(9) * t237;
t297 = t233 * t205 + (-t197 - t206) * t232 + t300;
t296 = t255 * t206 + (-t227 - t241) * t233 + t299;
t285 = Icges(2,4) * t289;
t275 = rSges(2,1) * t289 - rSges(2,2) * t287;
t274 = rSges(2,1) * t287 + rSges(2,2) * t289;
t273 = Icges(2,1) * t289 - t327;
t272 = Icges(2,1) * t287 + t285;
t271 = -Icges(2,2) * t287 + t285;
t270 = Icges(2,2) * t289 + t327;
t265 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t264 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t263 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t247 = t290 * rSges(3,3) + (rSges(3,1) * t293 + rSges(3,2) * t295) * t288;
t246 = Icges(3,5) * t290 + (Icges(3,1) * t293 + Icges(3,4) * t295) * t288;
t245 = Icges(3,6) * t290 + (Icges(3,4) * t293 + Icges(3,2) * t295) * t288;
t244 = Icges(3,3) * t290 + (Icges(3,5) * t293 + Icges(3,6) * t295) * t288;
t243 = V_base(5) * rSges(2,3) - t274 * V_base(6) + t314;
t242 = t275 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t231 = t274 * V_base(4) - t275 * V_base(5) + t310;
t228 = qJD(5) * t259 + t255;
t224 = t259 * rSges(4,1) - t258 * rSges(4,2) - rSges(4,3) * t321;
t223 = -rSges(5,1) * t321 - t259 * rSges(5,2) + t258 * rSges(5,3);
t216 = rSges(3,1) * t254 - rSges(3,2) * t253 + rSges(3,3) * t323;
t215 = rSges(3,1) * t252 - rSges(3,2) * t251 - rSges(3,3) * t322;
t214 = Icges(3,1) * t254 - Icges(3,4) * t253 + Icges(3,5) * t323;
t213 = Icges(3,1) * t252 - Icges(3,4) * t251 - Icges(3,5) * t322;
t212 = Icges(3,4) * t254 - Icges(3,2) * t253 + Icges(3,6) * t323;
t211 = Icges(3,4) * t252 - Icges(3,2) * t251 - Icges(3,6) * t322;
t210 = Icges(3,5) * t254 - Icges(3,6) * t253 + Icges(3,3) * t323;
t209 = Icges(3,5) * t252 - Icges(3,6) * t251 - Icges(3,3) * t322;
t199 = qJD(5) * t237 + t233;
t198 = qJD(5) * t235 + t232;
t192 = rSges(6,1) * t239 + rSges(6,2) * t238 + rSges(6,3) * t259;
t182 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t253;
t181 = rSges(4,1) * t235 - rSges(4,2) * t234 + rSges(4,3) * t251;
t180 = rSges(5,1) * t253 - rSges(5,2) * t237 + rSges(5,3) * t236;
t179 = rSges(5,1) * t251 - rSges(5,2) * t235 + rSges(5,3) * t234;
t165 = -t215 * t279 + t247 * t266 + t307;
t164 = t216 * t279 - t247 * t267 + t305;
t161 = rSges(6,1) * t204 + rSges(6,2) * t203 + rSges(6,3) * t237;
t159 = rSges(6,1) * t202 + rSges(6,2) * t201 + rSges(6,3) * t235;
t145 = t215 * t267 - t216 * t266 + t306;
t144 = -t181 * t255 + t224 * t232 + t304;
t143 = t182 * t255 - t224 * t233 + t302;
t142 = t181 * t233 - t182 * t232 + t303;
t141 = t223 * t232 + (-t179 - t196) * t255 + t301;
t140 = t180 * t255 + (-t223 - t227) * t233 + t299;
t139 = t179 * t233 + (-t180 - t197) * t232 + t300;
t138 = -t159 * t228 + t192 * t198 + t298;
t137 = t161 * t228 - t192 * t199 + t296;
t136 = t159 * t199 - t161 * t198 + t297;
t135 = qJD(6) * t237 + t198 * t316 - t228 * t318 + t298;
t134 = qJD(6) * t235 - t199 * t316 + t228 * t317 + t296;
t133 = qJD(6) * t259 - t198 * t317 + t199 * t318 + t297;
t1 = t279 * ((t209 * t266 + t210 * t267 + t244 * t279) * t290 + ((t212 * t295 + t214 * t293) * t267 + (t211 * t295 + t213 * t293) * t266 + (t245 * t295 + t246 * t293) * t279) * t288) / 0.2e1 + m(3) * (t145 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t267 * ((t210 * t323 - t253 * t212 + t254 * t214) * t267 + (t209 * t323 - t211 * t253 + t213 * t254) * t266 + (t244 * t323 - t245 * t253 + t246 * t254) * t279) / 0.2e1 + t266 * ((-t210 * t322 - t212 * t251 + t214 * t252) * t267 + (-t209 * t322 - t251 * t211 + t252 * t213) * t266 + (-t244 * t322 - t245 * t251 + t246 * t252) * t279) / 0.2e1 + m(2) * (t231 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + m(1) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + ((t201 * t344 + t202 * t343 + t235 * t345) * t228 + (t201 * t354 + t202 * t352 + t235 * t356) * t199 + (t355 * t201 + t353 * t202 + t357 * t235) * t198) * t198 / 0.2e1 + ((t203 * t344 + t204 * t343 + t237 * t345) * t228 + (t354 * t203 + t352 * t204 + t356 * t237) * t199 + (t203 * t355 + t204 * t353 + t237 * t357) * t198) * t199 / 0.2e1 + ((t344 * t238 + t343 * t239 + t345 * t259) * t228 + (t238 * t354 + t239 * t352 + t259 * t356) * t199 + (t238 * t355 + t239 * t353 + t259 * t357) * t198) * t228 / 0.2e1 + ((t234 * t342 + t235 * t341 + t251 * t340) * t255 + (t234 * t350 + t235 * t346 + t251 * t348) * t233 + (t351 * t234 + t347 * t235 + t349 * t251) * t232) * t232 / 0.2e1 + ((t236 * t342 + t237 * t341 + t253 * t340) * t255 + (t350 * t236 + t346 * t237 + t348 * t253) * t233 + (t236 * t351 + t237 * t347 + t253 * t349) * t232) * t233 / 0.2e1 + ((t342 * t258 + t341 * t259 - t340 * t321) * t255 + (t258 * t350 + t259 * t346 - t321 * t348) * t233 + (t258 * t351 + t259 * t347 - t321 * t349) * t232) * t255 / 0.2e1 + ((-t270 * t287 + t272 * t289 + Icges(1,4)) * V_base(5) + (-t287 * t271 + t289 * t273 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t289 * t270 + t287 * t272 + Icges(1,2)) * V_base(5) + (t271 * t289 + t273 * t287 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t287 + Icges(2,6) * t289 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t289 - Icges(2,6) * t287 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
