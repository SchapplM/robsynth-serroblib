% Calculate kinetic energy for
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:10:52
% EndTime: 2019-03-08 23:10:56
% DurationCPUTime: 4.49s
% Computational Cost: add. (3097->397), mult. (5456->572), div. (0->0), fcn. (6376->12), ass. (0->177)
t365 = Icges(5,1) + Icges(6,2);
t364 = Icges(6,1) + Icges(5,3);
t363 = -Icges(5,4) - Icges(6,6);
t362 = Icges(6,4) - Icges(5,5);
t361 = Icges(6,5) - Icges(5,6);
t360 = Icges(5,2) + Icges(6,3);
t298 = sin(pkin(11));
t300 = cos(pkin(11));
t307 = cos(qJ(2));
t301 = cos(pkin(6));
t304 = sin(qJ(2));
t335 = t301 * t304;
t263 = t298 * t307 + t300 * t335;
t299 = sin(pkin(6));
t333 = qJ(3) + qJ(4);
t324 = cos(t333);
t322 = t299 * t324;
t323 = sin(t333);
t233 = t263 * t323 + t300 * t322;
t321 = t299 * t323;
t234 = t263 * t324 - t300 * t321;
t334 = t301 * t307;
t262 = t298 * t304 - t300 * t334;
t358 = t360 * t233 + t363 * t234 + t361 * t262;
t265 = -t298 * t335 + t300 * t307;
t235 = t265 * t323 - t298 * t322;
t236 = t265 * t324 + t298 * t321;
t264 = t298 * t334 + t300 * t304;
t357 = t360 * t235 + t363 * t236 + t361 * t264;
t356 = t361 * t233 - t362 * t234 + t364 * t262;
t355 = t361 * t235 - t362 * t236 + t364 * t264;
t354 = t363 * t233 + t365 * t234 - t362 * t262;
t353 = t363 * t235 + t365 * t236 - t362 * t264;
t256 = -t301 * t324 + t304 * t321;
t257 = t301 * t323 + t304 * t322;
t337 = t299 * t307;
t352 = t360 * t256 + t363 * t257 - t361 * t337;
t351 = t363 * t256 + t365 * t257 + t362 * t337;
t350 = t361 * t256 - t362 * t257 - t364 * t337;
t345 = pkin(7) * t301;
t306 = cos(qJ(3));
t344 = pkin(3) * t306;
t342 = Icges(2,4) * t298;
t341 = t298 * t299;
t340 = t299 * t300;
t303 = sin(qJ(3));
t339 = t299 * t303;
t338 = t299 * t306;
t336 = t301 * t303;
t332 = qJD(2) * t299;
t331 = V_base(5) * qJ(1) + V_base(1);
t327 = qJD(1) + V_base(3);
t326 = t298 * t339;
t325 = t300 * t339;
t278 = t298 * t332 + V_base(4);
t290 = qJD(2) * t301 + V_base(6);
t241 = qJD(3) * t264 + t278;
t212 = qJD(4) * t264 + t241;
t277 = -t300 * t332 + V_base(5);
t240 = qJD(3) * t262 + t277;
t272 = pkin(1) * t298 - pkin(7) * t340;
t320 = -t272 * V_base(6) + V_base(5) * t345 + t331;
t273 = pkin(1) * t300 + pkin(7) * t341;
t319 = V_base(4) * t272 - t273 * V_base(5) + t327;
t211 = qJD(4) * t262 + t240;
t250 = (-qJD(3) - qJD(4)) * t337 + t290;
t318 = V_base(6) * t273 + V_base(2) + (-qJ(1) - t345) * V_base(4);
t229 = t263 * pkin(2) + t262 * pkin(8);
t271 = (pkin(2) * t304 - pkin(8) * t307) * t299;
t317 = -t229 * t290 + t277 * t271 + t320;
t230 = t265 * pkin(2) + t264 * pkin(8);
t316 = t278 * t229 - t230 * t277 + t319;
t315 = t290 * t230 - t271 * t278 + t318;
t181 = -pkin(3) * t325 + pkin(9) * t262 + t263 * t344;
t228 = pkin(3) * t336 + (-pkin(9) * t307 + t304 * t344) * t299;
t266 = -qJD(3) * t337 + t290;
t314 = -t181 * t266 + t240 * t228 + t317;
t182 = pkin(3) * t326 + pkin(9) * t264 + t265 * t344;
t313 = t241 * t181 - t182 * t240 + t316;
t312 = t266 * t182 - t228 * t241 + t315;
t226 = pkin(4) * t257 + qJ(5) * t256;
t311 = qJD(5) * t235 + t211 * t226 + t314;
t194 = pkin(4) * t234 + qJ(5) * t233;
t310 = qJD(5) * t256 + t212 * t194 + t313;
t195 = pkin(4) * t236 + qJ(5) * t235;
t309 = qJD(5) * t233 + t250 * t195 + t312;
t305 = cos(qJ(6));
t302 = sin(qJ(6));
t296 = Icges(2,4) * t300;
t288 = rSges(2,1) * t300 - rSges(2,2) * t298;
t287 = rSges(2,1) * t298 + rSges(2,2) * t300;
t286 = Icges(2,1) * t300 - t342;
t285 = Icges(2,1) * t298 + t296;
t284 = -Icges(2,2) * t298 + t296;
t283 = Icges(2,2) * t300 + t342;
t276 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t275 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t274 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t270 = t304 * t338 + t336;
t269 = t301 * t306 - t304 * t339;
t255 = rSges(3,3) * t301 + (rSges(3,1) * t304 + rSges(3,2) * t307) * t299;
t254 = Icges(3,5) * t301 + (Icges(3,1) * t304 + Icges(3,4) * t307) * t299;
t253 = Icges(3,6) * t301 + (Icges(3,4) * t304 + Icges(3,2) * t307) * t299;
t252 = Icges(3,3) * t301 + (Icges(3,5) * t304 + Icges(3,6) * t307) * t299;
t249 = V_base(5) * rSges(2,3) - t287 * V_base(6) + t331;
t248 = t288 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t246 = -pkin(5) * t337 + pkin(10) * t257;
t245 = t265 * t306 + t326;
t244 = -t265 * t303 + t298 * t338;
t243 = t263 * t306 - t325;
t242 = -t263 * t303 - t300 * t338;
t239 = t287 * V_base(4) - t288 * V_base(5) + t327;
t238 = t256 * t302 - t305 * t337;
t237 = t256 * t305 + t302 * t337;
t227 = rSges(4,1) * t270 + rSges(4,2) * t269 - rSges(4,3) * t337;
t225 = Icges(4,1) * t270 + Icges(4,4) * t269 - Icges(4,5) * t337;
t224 = Icges(4,4) * t270 + Icges(4,2) * t269 - Icges(4,6) * t337;
t223 = Icges(4,5) * t270 + Icges(4,6) * t269 - Icges(4,3) * t337;
t222 = rSges(3,1) * t265 - rSges(3,2) * t264 + rSges(3,3) * t341;
t221 = rSges(3,1) * t263 - rSges(3,2) * t262 - rSges(3,3) * t340;
t220 = Icges(3,1) * t265 - Icges(3,4) * t264 + Icges(3,5) * t341;
t219 = Icges(3,1) * t263 - Icges(3,4) * t262 - Icges(3,5) * t340;
t218 = Icges(3,4) * t265 - Icges(3,2) * t264 + Icges(3,6) * t341;
t217 = Icges(3,4) * t263 - Icges(3,2) * t262 - Icges(3,6) * t340;
t216 = Icges(3,5) * t265 - Icges(3,6) * t264 + Icges(3,3) * t341;
t215 = Icges(3,5) * t263 - Icges(3,6) * t262 - Icges(3,3) * t340;
t213 = qJD(6) * t257 + t250;
t209 = rSges(5,1) * t257 - rSges(5,2) * t256 - rSges(5,3) * t337;
t208 = -rSges(6,1) * t337 - rSges(6,2) * t257 + rSges(6,3) * t256;
t201 = pkin(5) * t264 + pkin(10) * t236;
t200 = pkin(5) * t262 + pkin(10) * t234;
t199 = t235 * t302 + t264 * t305;
t198 = t235 * t305 - t264 * t302;
t197 = t233 * t302 + t262 * t305;
t196 = t233 * t305 - t262 * t302;
t192 = rSges(4,1) * t245 + rSges(4,2) * t244 + rSges(4,3) * t264;
t191 = rSges(4,1) * t243 + rSges(4,2) * t242 + rSges(4,3) * t262;
t190 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t264;
t189 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t262;
t188 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t264;
t187 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t262;
t186 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t264;
t185 = Icges(4,5) * t243 + Icges(4,6) * t242 + Icges(4,3) * t262;
t184 = qJD(6) * t236 + t212;
t183 = qJD(6) * t234 + t211;
t180 = rSges(5,1) * t236 - rSges(5,2) * t235 + rSges(5,3) * t264;
t179 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t262;
t178 = rSges(6,1) * t264 - rSges(6,2) * t236 + rSges(6,3) * t235;
t177 = rSges(6,1) * t262 - rSges(6,2) * t234 + rSges(6,3) * t233;
t164 = rSges(7,1) * t238 + rSges(7,2) * t237 + rSges(7,3) * t257;
t163 = Icges(7,1) * t238 + Icges(7,4) * t237 + Icges(7,5) * t257;
t162 = Icges(7,4) * t238 + Icges(7,2) * t237 + Icges(7,6) * t257;
t161 = Icges(7,5) * t238 + Icges(7,6) * t237 + Icges(7,3) * t257;
t157 = -t221 * t290 + t255 * t277 + t320;
t156 = t222 * t290 - t255 * t278 + t318;
t153 = t221 * t278 - t222 * t277 + t319;
t152 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t236;
t151 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t234;
t150 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t236;
t149 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t234;
t148 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t236;
t147 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t234;
t146 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t236;
t145 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t234;
t144 = -t191 * t266 + t227 * t240 + t317;
t143 = t192 * t266 - t227 * t241 + t315;
t142 = t191 * t241 - t192 * t240 + t316;
t141 = -t179 * t250 + t209 * t211 + t314;
t140 = t180 * t250 - t209 * t212 + t312;
t139 = t179 * t212 - t180 * t211 + t313;
t138 = t208 * t211 + (-t177 - t194) * t250 + t311;
t137 = t178 * t250 + (-t208 - t226) * t212 + t309;
t136 = t177 * t212 + (-t178 - t195) * t211 + t310;
t135 = (-t194 - t200) * t250 - t151 * t213 + t164 * t183 + t211 * t246 + t311;
t134 = t152 * t213 - t164 * t184 + t201 * t250 + (-t226 - t246) * t212 + t309;
t133 = t151 * t184 - t152 * t183 + t200 * t212 + (-t195 - t201) * t211 + t310;
t1 = t278 * ((t216 * t341 - t218 * t264 + t220 * t265) * t278 + (t215 * t341 - t217 * t264 + t219 * t265) * t277 + (t252 * t341 - t253 * t264 + t254 * t265) * t290) / 0.2e1 + t277 * ((-t216 * t340 - t218 * t262 + t220 * t263) * t278 + (-t215 * t340 - t262 * t217 + t263 * t219) * t277 + (-t252 * t340 - t253 * t262 + t254 * t263) * t290) / 0.2e1 + t266 * ((-t186 * t337 + t188 * t269 + t190 * t270) * t241 + (-t185 * t337 + t187 * t269 + t189 * t270) * t240 + (-t223 * t337 + t269 * t224 + t270 * t225) * t266) / 0.2e1 + m(3) * (t153 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + t290 * ((t215 * t277 + t216 * t278 + t252 * t290) * t301 + ((t218 * t307 + t220 * t304) * t278 + (t217 * t307 + t219 * t304) * t277 + (t253 * t307 + t254 * t304) * t290) * t299) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t183 * ((t146 * t234 + t148 * t196 + t150 * t197) * t184 + (t234 * t145 + t196 * t147 + t197 * t149) * t183 + (t161 * t234 + t162 * t196 + t163 * t197) * t213) / 0.2e1 + t184 * ((t236 * t146 + t198 * t148 + t199 * t150) * t184 + (t145 * t236 + t147 * t198 + t149 * t199) * t183 + (t161 * t236 + t162 * t198 + t163 * t199) * t213) / 0.2e1 + m(2) * (t239 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t213 * ((t146 * t257 + t148 * t237 + t150 * t238) * t184 + (t145 * t257 + t147 * t237 + t149 * t238) * t183 + (t161 * t257 + t237 * t162 + t238 * t163) * t213) / 0.2e1 + t240 * ((t186 * t262 + t188 * t242 + t190 * t243) * t241 + (t262 * t185 + t242 * t187 + t189 * t243) * t240 + (t223 * t262 + t224 * t242 + t225 * t243) * t266) / 0.2e1 + t241 * ((t264 * t186 + t244 * t188 + t245 * t190) * t241 + (t185 * t264 + t187 * t244 + t189 * t245) * t240 + (t223 * t264 + t224 * t244 + t225 * t245) * t266) / 0.2e1 + m(1) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + ((t233 * t352 + t234 * t351 + t262 * t350) * t250 + (t233 * t357 + t234 * t353 + t262 * t355) * t212 + (t358 * t233 + t354 * t234 + t356 * t262) * t211) * t211 / 0.2e1 + ((t235 * t352 + t236 * t351 + t264 * t350) * t250 + (t357 * t235 + t353 * t236 + t355 * t264) * t212 + (t235 * t358 + t236 * t354 + t264 * t356) * t211) * t212 / 0.2e1 + ((t352 * t256 + t351 * t257 - t350 * t337) * t250 + (t256 * t357 + t257 * t353 - t337 * t355) * t212 + (t256 * t358 + t257 * t354 - t337 * t356) * t211) * t250 / 0.2e1 + ((-t283 * t298 + t285 * t300 + Icges(1,4)) * V_base(5) + (-t298 * t284 + t300 * t286 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t283 + t298 * t285 + Icges(1,2)) * V_base(5) + (t284 * t300 + t286 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t300 - Icges(2,6) * t298 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t298 + Icges(2,6) * t300 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
