% Calculate kinetic energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:33
% EndTime: 2019-03-08 21:22:37
% DurationCPUTime: 4.31s
% Computational Cost: add. (3452->386), mult. (5968->540), div. (0->0), fcn. (7093->12), ass. (0->174)
t384 = Icges(6,1) + Icges(7,1);
t383 = Icges(6,4) + Icges(7,4);
t382 = Icges(6,5) + Icges(7,5);
t381 = Icges(6,2) + Icges(7,2);
t380 = Icges(6,6) + Icges(7,6);
t379 = Icges(4,3) + Icges(5,3);
t378 = Icges(6,3) + Icges(7,3);
t377 = rSges(7,3) + qJ(6);
t301 = sin(pkin(10));
t303 = cos(pkin(10));
t312 = cos(qJ(2));
t304 = cos(pkin(6));
t309 = sin(qJ(2));
t341 = t304 * t309;
t265 = t301 * t312 + t303 * t341;
t335 = qJ(3) + pkin(11);
t298 = sin(t335);
t326 = cos(t335);
t302 = sin(pkin(6));
t347 = t302 * t303;
t239 = t265 * t326 - t298 * t347;
t340 = t304 * t312;
t264 = t301 * t309 - t303 * t340;
t307 = sin(qJ(5));
t310 = cos(qJ(5));
t211 = -t239 * t307 + t264 * t310;
t350 = t264 * t307;
t212 = t239 * t310 + t350;
t325 = t302 * t326;
t238 = t265 * t298 + t303 * t325;
t376 = t211 * t380 + t212 * t382 + t238 * t378;
t267 = -t301 * t341 + t303 * t312;
t348 = t301 * t302;
t241 = t267 * t326 + t298 * t348;
t266 = t301 * t340 + t303 * t309;
t213 = -t241 * t307 + t266 * t310;
t349 = t266 * t307;
t214 = t241 * t310 + t349;
t240 = t267 * t298 - t301 * t325;
t375 = t213 * t380 + t214 * t382 + t240 * t378;
t374 = t211 * t381 + t212 * t383 + t238 * t380;
t373 = t213 * t381 + t214 * t383 + t240 * t380;
t372 = t211 * t383 + t212 * t384 + t238 * t382;
t371 = t213 * t383 + t214 * t384 + t240 * t382;
t258 = t298 * t304 + t309 * t325;
t343 = t302 * t312;
t242 = -t258 * t307 - t310 * t343;
t327 = t307 * t343;
t243 = t258 * t310 - t327;
t345 = t302 * t309;
t257 = t298 * t345 - t304 * t326;
t370 = t242 * t380 + t243 * t382 + t257 * t378;
t369 = t242 * t381 + t243 * t383 + t257 * t380;
t368 = t242 * t383 + t243 * t384 + t257 * t382;
t308 = sin(qJ(3));
t311 = cos(qJ(3));
t344 = t302 * t311;
t247 = -t265 * t308 - t303 * t344;
t346 = t302 * t308;
t328 = t303 * t346;
t248 = t265 * t311 - t328;
t367 = Icges(4,5) * t248 + Icges(5,5) * t239 + Icges(4,6) * t247 - Icges(5,6) * t238 + t264 * t379;
t249 = -t267 * t308 + t301 * t344;
t329 = t301 * t346;
t250 = t267 * t311 + t329;
t366 = Icges(4,5) * t250 + Icges(5,5) * t241 + Icges(4,6) * t249 - Icges(5,6) * t240 + t266 * t379;
t271 = t304 * t311 - t308 * t345;
t342 = t304 * t308;
t272 = t309 * t344 + t342;
t365 = Icges(4,5) * t272 + Icges(5,5) * t258 + Icges(4,6) * t271 - Icges(5,6) * t257 - t343 * t379;
t356 = pkin(7) * t304;
t355 = pkin(3) * t311;
t354 = pkin(5) * t310;
t351 = Icges(2,4) * t301;
t339 = rSges(7,1) * t212 + rSges(7,2) * t211 + pkin(5) * t350 + t238 * t377 + t239 * t354;
t338 = rSges(7,1) * t214 + rSges(7,2) * t213 + pkin(5) * t349 + t240 * t377 + t241 * t354;
t337 = rSges(7,1) * t243 + rSges(7,2) * t242 - pkin(5) * t327 + t257 * t377 + t258 * t354;
t336 = qJD(2) * t302;
t334 = V_base(5) * qJ(1) + V_base(1);
t330 = qJD(1) + V_base(3);
t281 = t301 * t336 + V_base(4);
t292 = qJD(2) * t304 + V_base(6);
t246 = qJD(3) * t266 + t281;
t280 = -t303 * t336 + V_base(5);
t245 = qJD(3) * t264 + t280;
t268 = -qJD(3) * t343 + t292;
t274 = pkin(1) * t301 - pkin(7) * t347;
t324 = -t274 * V_base(6) + t356 * V_base(5) + t334;
t275 = pkin(1) * t303 + pkin(7) * t348;
t323 = t274 * V_base(4) - t275 * V_base(5) + t330;
t322 = V_base(6) * t275 + V_base(2) + (-qJ(1) - t356) * V_base(4);
t236 = pkin(2) * t265 + pkin(8) * t264;
t273 = (pkin(2) * t309 - pkin(8) * t312) * t302;
t321 = -t236 * t292 + t273 * t280 + t324;
t237 = pkin(2) * t267 + pkin(8) * t266;
t320 = t236 * t281 - t237 * t280 + t323;
t319 = t237 * t292 - t273 * t281 + t322;
t233 = pkin(3) * t342 + (-qJ(4) * t312 + t309 * t355) * t302;
t318 = qJD(4) * t266 + t233 * t245 + t321;
t195 = pkin(3) * t329 + qJ(4) * t266 + t267 * t355;
t317 = qJD(4) * t264 + t195 * t268 + t319;
t194 = -pkin(3) * t328 + qJ(4) * t264 + t265 * t355;
t316 = -qJD(4) * t343 + t194 * t246 + t320;
t207 = pkin(4) * t239 + pkin(9) * t238;
t227 = pkin(4) * t258 + pkin(9) * t257;
t315 = t245 * t227 + (-t194 - t207) * t268 + t318;
t208 = pkin(4) * t241 + pkin(9) * t240;
t314 = t268 * t208 + (-t227 - t233) * t246 + t317;
t313 = t246 * t207 + (-t195 - t208) * t245 + t316;
t299 = Icges(2,4) * t303;
t289 = rSges(2,1) * t303 - rSges(2,2) * t301;
t288 = rSges(2,1) * t301 + rSges(2,2) * t303;
t287 = Icges(2,1) * t303 - t351;
t286 = Icges(2,1) * t301 + t299;
t285 = -Icges(2,2) * t301 + t299;
t284 = Icges(2,2) * t303 + t351;
t279 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t278 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t277 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t259 = t304 * rSges(3,3) + (rSges(3,1) * t309 + rSges(3,2) * t312) * t302;
t256 = Icges(3,5) * t304 + (Icges(3,1) * t309 + Icges(3,4) * t312) * t302;
t255 = Icges(3,6) * t304 + (Icges(3,4) * t309 + Icges(3,2) * t312) * t302;
t254 = Icges(3,3) * t304 + (Icges(3,5) * t309 + Icges(3,6) * t312) * t302;
t253 = V_base(5) * rSges(2,3) - t288 * V_base(6) + t334;
t252 = t289 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t244 = t288 * V_base(4) - t289 * V_base(5) + t330;
t235 = qJD(5) * t257 + t268;
t234 = rSges(4,1) * t272 + rSges(4,2) * t271 - rSges(4,3) * t343;
t232 = Icges(4,1) * t272 + Icges(4,4) * t271 - Icges(4,5) * t343;
t231 = Icges(4,4) * t272 + Icges(4,2) * t271 - Icges(4,6) * t343;
t229 = rSges(3,1) * t267 - rSges(3,2) * t266 + rSges(3,3) * t348;
t228 = rSges(3,1) * t265 - rSges(3,2) * t264 - rSges(3,3) * t347;
t226 = Icges(3,1) * t267 - Icges(3,4) * t266 + Icges(3,5) * t348;
t225 = Icges(3,1) * t265 - Icges(3,4) * t264 - Icges(3,5) * t347;
t224 = Icges(3,4) * t267 - Icges(3,2) * t266 + Icges(3,6) * t348;
t223 = Icges(3,4) * t265 - Icges(3,2) * t264 - Icges(3,6) * t347;
t222 = Icges(3,5) * t267 - Icges(3,6) * t266 + Icges(3,3) * t348;
t221 = Icges(3,5) * t265 - Icges(3,6) * t264 - Icges(3,3) * t347;
t218 = rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t343;
t217 = Icges(5,1) * t258 - Icges(5,4) * t257 - Icges(5,5) * t343;
t216 = Icges(5,4) * t258 - Icges(5,2) * t257 - Icges(5,6) * t343;
t210 = qJD(5) * t240 + t246;
t209 = qJD(5) * t238 + t245;
t204 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t266;
t203 = rSges(4,1) * t248 + rSges(4,2) * t247 + rSges(4,3) * t264;
t202 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t266;
t201 = Icges(4,1) * t248 + Icges(4,4) * t247 + Icges(4,5) * t264;
t200 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t266;
t199 = Icges(4,4) * t248 + Icges(4,2) * t247 + Icges(4,6) * t264;
t193 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t266;
t192 = rSges(5,1) * t239 - rSges(5,2) * t238 + rSges(5,3) * t264;
t191 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t266;
t190 = Icges(5,1) * t239 - Icges(5,4) * t238 + Icges(5,5) * t264;
t189 = Icges(5,4) * t241 - Icges(5,2) * t240 + Icges(5,6) * t266;
t188 = Icges(5,4) * t239 - Icges(5,2) * t238 + Icges(5,6) * t264;
t185 = rSges(6,1) * t243 + rSges(6,2) * t242 + rSges(6,3) * t257;
t174 = -t228 * t292 + t259 * t280 + t324;
t173 = t229 * t292 - t259 * t281 + t322;
t171 = t228 * t281 - t229 * t280 + t323;
t170 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t240;
t168 = rSges(6,1) * t212 + rSges(6,2) * t211 + rSges(6,3) * t238;
t152 = -t203 * t268 + t234 * t245 + t321;
t151 = t204 * t268 - t234 * t246 + t319;
t150 = t203 * t246 - t204 * t245 + t320;
t149 = t218 * t245 + (-t192 - t194) * t268 + t318;
t148 = t193 * t268 + (-t218 - t233) * t246 + t317;
t147 = t246 * t192 + (-t193 - t195) * t245 + t316;
t146 = -t168 * t235 + t185 * t209 + t315;
t145 = t170 * t235 - t185 * t210 + t314;
t144 = t168 * t210 - t170 * t209 + t313;
t143 = qJD(6) * t240 + t209 * t337 - t235 * t339 + t315;
t142 = qJD(6) * t238 - t210 * t337 + t235 * t338 + t314;
t141 = qJD(6) * t257 - t209 * t338 + t210 * t339 + t313;
t1 = m(4) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(6) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + t292 * ((t221 * t280 + t222 * t281 + t254 * t292) * t304 + ((t224 * t312 + t226 * t309) * t281 + (t223 * t312 + t225 * t309) * t280 + (t255 * t312 + t256 * t309) * t292) * t302) / 0.2e1 + m(7) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(3) * (t171 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(2) * (t244 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + t280 * ((-t222 * t347 - t224 * t264 + t226 * t265) * t281 + (-t221 * t347 - t223 * t264 + t225 * t265) * t280 + (-t254 * t347 - t255 * t264 + t256 * t265) * t292) / 0.2e1 + t281 * ((t222 * t348 - t224 * t266 + t226 * t267) * t281 + (t221 * t348 - t223 * t266 + t225 * t267) * t280 + (t254 * t348 - t255 * t266 + t256 * t267) * t292) / 0.2e1 + m(1) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + ((t211 * t369 + t212 * t368 + t238 * t370) * t235 + (t211 * t373 + t212 * t371 + t238 * t375) * t210 + (t374 * t211 + t372 * t212 + t376 * t238) * t209) * t209 / 0.2e1 + ((t213 * t369 + t214 * t368 + t240 * t370) * t235 + (t373 * t213 + t371 * t214 + t375 * t240) * t210 + (t213 * t374 + t214 * t372 + t240 * t376) * t209) * t210 / 0.2e1 + ((t369 * t242 + t368 * t243 + t370 * t257) * t235 + (t242 * t373 + t243 * t371 + t257 * t375) * t210 + (t242 * t374 + t243 * t372 + t257 * t376) * t209) * t235 / 0.2e1 + ((-t216 * t238 + t217 * t239 + t231 * t247 + t232 * t248 + t264 * t365) * t268 + (-t189 * t238 + t191 * t239 + t200 * t247 + t202 * t248 + t264 * t366) * t246 + (-t188 * t238 + t190 * t239 + t199 * t247 + t201 * t248 + t367 * t264) * t245) * t245 / 0.2e1 + ((-t216 * t240 + t217 * t241 + t231 * t249 + t232 * t250 + t266 * t365) * t268 + (-t189 * t240 + t191 * t241 + t200 * t249 + t202 * t250 + t366 * t266) * t246 + (-t188 * t240 + t190 * t241 + t199 * t249 + t201 * t250 + t266 * t367) * t245) * t246 / 0.2e1 + ((-t257 * t216 + t258 * t217 + t271 * t231 + t272 * t232 - t365 * t343) * t268 + (-t257 * t189 + t258 * t191 + t271 * t200 + t272 * t202 - t343 * t366) * t246 + (-t257 * t188 + t258 * t190 + t271 * t199 + t272 * t201 - t343 * t367) * t245) * t268 / 0.2e1 + ((-t284 * t301 + t286 * t303 + Icges(1,4)) * V_base(5) + (-t285 * t301 + t287 * t303 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t284 * t303 + t286 * t301 + Icges(1,2)) * V_base(5) + (t285 * t303 + t287 * t301 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t303 - Icges(2,6) * t301 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t301 + Icges(2,6) * t303 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
