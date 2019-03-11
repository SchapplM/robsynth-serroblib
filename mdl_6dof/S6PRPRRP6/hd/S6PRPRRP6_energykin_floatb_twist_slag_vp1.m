% Calculate kinetic energy for
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:07
% EndTime: 2019-03-08 20:18:11
% DurationCPUTime: 4.01s
% Computational Cost: add. (2272->340), mult. (5309->474), div. (0->0), fcn. (6265->10), ass. (0->156)
t380 = Icges(3,1) + Icges(4,2);
t379 = Icges(6,1) + Icges(7,1);
t378 = Icges(3,4) + Icges(4,6);
t377 = -Icges(6,4) + Icges(7,5);
t376 = Icges(7,4) + Icges(6,5);
t375 = Icges(3,5) - Icges(4,4);
t374 = Icges(3,2) + Icges(4,3);
t373 = Icges(6,2) + Icges(7,3);
t372 = Icges(7,2) + Icges(6,3);
t371 = Icges(3,6) - Icges(4,5);
t370 = Icges(6,6) - Icges(7,6);
t369 = Icges(3,3) + Icges(4,1);
t368 = rSges(7,1) + pkin(5);
t367 = rSges(7,3) + qJ(6);
t300 = sin(pkin(10));
t302 = cos(pkin(10));
t306 = sin(qJ(2));
t303 = cos(pkin(6));
t307 = cos(qJ(2));
t330 = t303 * t307;
t266 = t300 * t330 + t302 * t306;
t305 = sin(qJ(4));
t301 = sin(pkin(6));
t340 = cos(qJ(4));
t320 = t301 * t340;
t238 = t266 * t305 + t300 * t320;
t331 = t303 * t306;
t267 = -t300 * t331 + t302 * t307;
t304 = sin(qJ(5));
t339 = cos(qJ(5));
t202 = t238 * t304 - t267 * t339;
t203 = t238 * t339 + t267 * t304;
t334 = t301 * t305;
t237 = -t266 * t340 + t300 * t334;
t366 = t373 * t202 + t377 * t203 - t370 * t237;
t264 = t300 * t306 - t302 * t330;
t240 = t264 * t305 - t302 * t320;
t265 = t300 * t307 + t302 * t331;
t204 = t240 * t304 - t265 * t339;
t205 = t240 * t339 + t265 * t304;
t239 = t264 * t340 + t302 * t334;
t365 = t373 * t204 + t377 * t205 + t370 * t239;
t364 = -t370 * t202 + t376 * t203 + t372 * t237;
t363 = -t370 * t204 + t376 * t205 - t372 * t239;
t362 = t377 * t202 + t379 * t203 + t376 * t237;
t361 = t377 * t204 + t379 * t205 - t376 * t239;
t332 = t301 * t307;
t272 = t303 * t340 - t305 * t332;
t333 = t301 * t306;
t241 = t272 * t304 - t333 * t339;
t242 = t272 * t339 + t304 * t333;
t271 = t303 * t305 + t307 * t320;
t360 = t373 * t241 + t377 * t242 - t370 * t271;
t359 = -t370 * t241 + t376 * t242 + t372 * t271;
t358 = t377 * t241 + t379 * t242 + t376 * t271;
t336 = t300 * t301;
t357 = t374 * t266 - t378 * t267 - t371 * t336;
t335 = t301 * t302;
t356 = t374 * t264 - t378 * t265 + t371 * t335;
t355 = -t378 * t266 + t380 * t267 + t375 * t336;
t354 = -t378 * t264 + t380 * t265 - t375 * t335;
t353 = -t371 * t266 + t375 * t267 + t369 * t336;
t352 = -t371 * t264 + t375 * t265 - t369 * t335;
t351 = t369 * t303 + (t375 * t306 + t371 * t307) * t301;
t350 = t371 * t303 + (t378 * t306 + t374 * t307) * t301;
t349 = t375 * t303 + (t380 * t306 + t378 * t307) * t301;
t338 = pkin(7) * t303;
t337 = Icges(2,4) * t300;
t329 = rSges(7,2) * t237 + t367 * t202 + t368 * t203;
t328 = -rSges(7,2) * t239 + t367 * t204 + t368 * t205;
t327 = rSges(7,2) * t271 + t367 * t241 + t368 * t242;
t326 = qJD(2) * t301;
t325 = V_base(5) * qJ(1) + V_base(1);
t321 = qJD(1) + V_base(3);
t281 = t300 * t326 + V_base(4);
t293 = qJD(2) * t303 + V_base(6);
t236 = qJD(4) * t267 + t281;
t268 = qJD(4) * t333 + t293;
t280 = -t302 * t326 + V_base(5);
t235 = qJD(4) * t265 + t280;
t274 = pkin(1) * t300 - pkin(7) * t335;
t319 = -t274 * V_base(6) + V_base(5) * t338 + t325;
t275 = pkin(1) * t302 + pkin(7) * t336;
t318 = V_base(4) * t274 - V_base(5) * t275 + t321;
t317 = V_base(6) * t275 + V_base(2) + (-qJ(1) - t338) * V_base(4);
t273 = (pkin(2) * t306 - qJ(3) * t307) * t301;
t316 = qJD(3) * t266 + t280 * t273 + t319;
t230 = pkin(2) * t267 + qJ(3) * t266;
t315 = qJD(3) * t264 + t293 * t230 + t317;
t229 = pkin(2) * t265 + qJ(3) * t264;
t314 = -qJD(3) * t332 + t281 * t229 + t318;
t245 = -pkin(3) * t335 + pkin(8) * t265;
t276 = pkin(3) * t303 + pkin(8) * t333;
t313 = t280 * t276 + (-t229 - t245) * t293 + t316;
t244 = pkin(3) * t336 + pkin(8) * t267;
t312 = t293 * t244 + (-t273 - t276) * t281 + t315;
t311 = t281 * t245 + (-t230 - t244) * t280 + t314;
t198 = pkin(4) * t240 - pkin(9) * t239;
t231 = pkin(4) * t272 + pkin(9) * t271;
t310 = -t198 * t268 + t235 * t231 + t313;
t197 = pkin(4) * t238 + pkin(9) * t237;
t309 = t268 * t197 - t231 * t236 + t312;
t308 = -t235 * t197 + t236 * t198 + t311;
t298 = Icges(2,4) * t302;
t289 = rSges(2,1) * t302 - rSges(2,2) * t300;
t288 = rSges(2,1) * t300 + rSges(2,2) * t302;
t287 = Icges(2,1) * t302 - t337;
t286 = Icges(2,1) * t300 + t298;
t285 = -Icges(2,2) * t300 + t298;
t284 = Icges(2,2) * t302 + t337;
t279 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t278 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t277 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t256 = t303 * rSges(4,1) + (-rSges(4,2) * t306 - rSges(4,3) * t307) * t301;
t255 = t303 * rSges(3,3) + (rSges(3,1) * t306 + rSges(3,2) * t307) * t301;
t247 = V_base(5) * rSges(2,3) - t288 * V_base(6) + t325;
t246 = t289 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t234 = t288 * V_base(4) - t289 * V_base(5) + t321;
t233 = qJD(5) * t271 + t268;
t227 = rSges(5,1) * t272 - rSges(5,2) * t271 + rSges(5,3) * t333;
t226 = Icges(5,1) * t272 - Icges(5,4) * t271 + Icges(5,5) * t333;
t225 = Icges(5,4) * t272 - Icges(5,2) * t271 + Icges(5,6) * t333;
t224 = Icges(5,5) * t272 - Icges(5,6) * t271 + Icges(5,3) * t333;
t223 = rSges(3,1) * t267 - rSges(3,2) * t266 + rSges(3,3) * t336;
t222 = rSges(3,1) * t265 - rSges(3,2) * t264 - rSges(3,3) * t335;
t221 = -rSges(4,1) * t335 - rSges(4,2) * t265 + rSges(4,3) * t264;
t220 = rSges(4,1) * t336 - rSges(4,2) * t267 + rSges(4,3) * t266;
t201 = qJD(5) * t237 + t236;
t200 = -qJD(5) * t239 + t235;
t195 = rSges(6,1) * t242 - rSges(6,2) * t241 + rSges(6,3) * t271;
t186 = rSges(5,1) * t240 + rSges(5,2) * t239 + rSges(5,3) * t265;
t185 = rSges(5,1) * t238 - rSges(5,2) * t237 + rSges(5,3) * t267;
t184 = Icges(5,1) * t240 + Icges(5,4) * t239 + Icges(5,5) * t265;
t183 = Icges(5,1) * t238 - Icges(5,4) * t237 + Icges(5,5) * t267;
t182 = Icges(5,4) * t240 + Icges(5,2) * t239 + Icges(5,6) * t265;
t181 = Icges(5,4) * t238 - Icges(5,2) * t237 + Icges(5,6) * t267;
t180 = Icges(5,5) * t240 + Icges(5,6) * t239 + Icges(5,3) * t265;
t179 = Icges(5,5) * t238 - Icges(5,6) * t237 + Icges(5,3) * t267;
t175 = -t222 * t293 + t255 * t280 + t319;
t174 = t223 * t293 - t255 * t281 + t317;
t173 = rSges(6,1) * t205 - rSges(6,2) * t204 - rSges(6,3) * t239;
t171 = rSges(6,1) * t203 - rSges(6,2) * t202 + rSges(6,3) * t237;
t157 = t222 * t281 - t223 * t280 + t318;
t156 = t256 * t280 + (-t221 - t229) * t293 + t316;
t155 = t220 * t293 + (-t256 - t273) * t281 + t315;
t154 = t281 * t221 + (-t220 - t230) * t280 + t314;
t153 = -t186 * t268 + t227 * t235 + t313;
t152 = t185 * t268 - t227 * t236 + t312;
t151 = -t235 * t185 + t236 * t186 + t311;
t150 = -t173 * t233 + t195 * t200 + t310;
t149 = t171 * t233 - t195 * t201 + t309;
t148 = -t200 * t171 + t201 * t173 + t308;
t147 = qJD(6) * t202 + t200 * t327 - t233 * t328 + t310;
t146 = qJD(6) * t204 - t201 * t327 + t233 * t329 + t309;
t145 = qJD(6) * t241 - t200 * t329 + t201 * t328 + t308;
t1 = t268 * ((t179 * t333 - t181 * t271 + t183 * t272) * t236 + (t180 * t333 - t182 * t271 + t184 * t272) * t235 + (t224 * t333 - t271 * t225 + t272 * t226) * t268) / 0.2e1 + m(7) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(6) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(5) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(4) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(3) * (t157 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(2) * (t234 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + t235 * ((t179 * t265 + t181 * t239 + t183 * t240) * t236 + (t265 * t180 + t239 * t182 + t240 * t184) * t235 + (t224 * t265 + t225 * t239 + t226 * t240) * t268) / 0.2e1 + t236 * ((t267 * t179 - t237 * t181 + t238 * t183) * t236 + (t180 * t267 - t182 * t237 + t184 * t238) * t235 + (t224 * t267 - t225 * t237 + t226 * t238) * t268) / 0.2e1 + m(1) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + ((t360 * t204 + t358 * t205 - t359 * t239) * t233 + (t366 * t204 + t362 * t205 - t364 * t239) * t201 + (t365 * t204 + t361 * t205 - t363 * t239) * t200) * t200 / 0.2e1 + ((t360 * t202 + t358 * t203 + t359 * t237) * t233 + (t366 * t202 + t362 * t203 + t364 * t237) * t201 + (t365 * t202 + t361 * t203 + t363 * t237) * t200) * t201 / 0.2e1 + ((t360 * t241 + t358 * t242 + t359 * t271) * t233 + (t366 * t241 + t362 * t242 + t364 * t271) * t201 + (t365 * t241 + t361 * t242 + t363 * t271) * t200) * t233 / 0.2e1 + ((-t350 * t264 + t349 * t265 - t351 * t335) * t293 + (t357 * t264 + t355 * t265 - t353 * t335) * t281 + (t356 * t264 + t354 * t265 - t352 * t335) * t280) * t280 / 0.2e1 + ((-t350 * t266 + t349 * t267 + t351 * t336) * t293 + (t357 * t266 + t355 * t267 + t353 * t336) * t281 + (t356 * t266 + t354 * t267 + t352 * t336) * t280) * t281 / 0.2e1 + ((t352 * t280 + t353 * t281 + t351 * t293) * t303 + ((t349 * t306 + t350 * t307) * t293 + (t355 * t306 - t357 * t307) * t281 + (t354 * t306 - t356 * t307) * t280) * t301) * t293 / 0.2e1 + ((-t284 * t300 + t286 * t302 + Icges(1,4)) * V_base(5) + (-t300 * t285 + t302 * t287 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t302 * t284 + t300 * t286 + Icges(1,2)) * V_base(5) + (t285 * t302 + t287 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t302 - Icges(2,6) * t300 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t300 + Icges(2,6) * t302 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
