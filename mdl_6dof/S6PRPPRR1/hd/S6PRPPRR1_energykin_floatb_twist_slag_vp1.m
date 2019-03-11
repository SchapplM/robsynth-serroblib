% Calculate kinetic energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:50
% EndTime: 2019-03-08 19:13:56
% DurationCPUTime: 5.86s
% Computational Cost: add. (3727->431), mult. (7711->604), div. (0->0), fcn. (9550->14), ass. (0->189)
t393 = -Icges(4,2) - Icges(5,3);
t392 = Icges(4,3) + Icges(3,3);
t332 = sin(qJ(2));
t375 = sin(pkin(11));
t376 = cos(pkin(11));
t380 = cos(qJ(2));
t293 = -t332 * t376 - t375 * t380;
t326 = sin(pkin(10));
t329 = cos(pkin(10));
t377 = cos(pkin(6));
t348 = t377 * t375;
t349 = t377 * t376;
t339 = -t332 * t348 + t349 * t380;
t263 = t326 * t293 + t329 * t339;
t285 = t332 * t349 + t348 * t380;
t294 = -t332 * t375 + t380 * t376;
t264 = t285 * t329 + t294 * t326;
t351 = t377 * t380;
t287 = -t326 * t332 + t329 * t351;
t355 = t332 * t377;
t288 = t326 * t380 + t329 * t355;
t327 = sin(pkin(6));
t372 = t327 * t329;
t389 = Icges(3,5) * t288 + Icges(4,5) * t264 + Icges(3,6) * t287 + Icges(4,6) * t263 - t392 * t372;
t265 = t293 * t329 - t326 * t339;
t266 = -t285 * t326 + t294 * t329;
t289 = -t326 * t351 - t329 * t332;
t290 = -t326 * t355 + t329 * t380;
t373 = t326 * t327;
t388 = Icges(3,5) * t290 + Icges(4,5) * t266 + Icges(3,6) * t289 + Icges(4,6) * t265 + t392 * t373;
t325 = sin(pkin(12));
t328 = cos(pkin(12));
t237 = -t264 * t325 - t328 * t372;
t358 = t325 * t372;
t238 = t264 * t328 - t358;
t387 = -Icges(4,4) * t264 + Icges(5,5) * t238 + Icges(4,6) * t372 + Icges(5,6) * t237 + t393 * t263;
t239 = -t266 * t325 + t328 * t373;
t359 = t325 * t373;
t240 = t266 * t328 + t359;
t386 = -Icges(4,4) * t266 + Icges(5,5) * t240 - Icges(4,6) * t373 + Icges(5,6) * t239 + t393 * t265;
t283 = t294 * t327;
t284 = t293 * t327;
t385 = -Icges(4,5) * t284 + Icges(4,6) * t283 + (Icges(3,5) * t332 + Icges(3,6) * t380) * t327 + t392 * t377;
t270 = t284 * t325 + t328 * t377;
t353 = t377 * t325;
t271 = -t284 * t328 + t353;
t384 = Icges(4,4) * t284 + Icges(5,5) * t271 - Icges(4,6) * t377 + Icges(5,6) * t270 + t393 * t283;
t379 = pkin(2) * t380;
t378 = pkin(4) * t328;
t374 = Icges(2,4) * t326;
t225 = pkin(3) * t264 - qJ(4) * t263;
t356 = pkin(2) * t355 - qJ(3) * t327;
t256 = t326 * t379 + t329 * t356;
t370 = -t225 - t256;
t226 = pkin(3) * t266 - qJ(4) * t265;
t257 = -t326 * t356 + t329 * t379;
t369 = -t226 - t257;
t258 = -t284 * pkin(3) - t283 * qJ(4);
t295 = t327 * t332 * pkin(2) + qJ(3) * t377;
t368 = -t258 - t295;
t367 = qJD(2) * t327;
t366 = qJD(3) * t327;
t365 = pkin(12) + qJ(5);
t364 = V_base(5) * qJ(1) + V_base(1);
t360 = qJD(1) + V_base(3);
t357 = t377 * pkin(7);
t303 = t326 * t367 + V_base(4);
t314 = qJD(2) * t377 + V_base(6);
t352 = cos(t365);
t242 = -qJD(5) * t265 + t303;
t272 = -qJD(5) * t283 + t314;
t347 = t327 * t352;
t302 = -t329 * t367 + V_base(5);
t241 = -qJD(5) * t263 + t302;
t296 = pkin(1) * t326 - pkin(7) * t372;
t346 = -t296 * V_base(6) + V_base(5) * t357 + t364;
t297 = pkin(1) * t329 + pkin(7) * t373;
t345 = V_base(4) * t296 - t297 * V_base(5) + t360;
t344 = t302 * t295 + t326 * t366 + t346;
t343 = qJD(3) * t377 + t303 * t256 + t345;
t342 = V_base(6) * t297 + V_base(2) + (-t357 - qJ(1)) * V_base(4);
t341 = -qJD(4) * t265 + t302 * t258 + t344;
t340 = -qJD(4) * t283 + t303 * t225 + t343;
t338 = t314 * t257 - t329 * t366 + t342;
t337 = -qJD(4) * t263 + t314 * t226 + t338;
t172 = -pkin(4) * t358 - pkin(8) * t263 + t264 * t378;
t206 = pkin(4) * t353 - pkin(8) * t283 - t284 * t378;
t336 = t302 * t206 + (-t172 + t370) * t314 + t341;
t173 = pkin(4) * t359 - pkin(8) * t265 + t266 * t378;
t335 = t303 * t172 + (-t173 + t369) * t302 + t340;
t334 = t314 * t173 + (-t206 + t368) * t303 + t337;
t333 = cos(qJ(6));
t331 = sin(qJ(6));
t323 = Icges(2,4) * t329;
t322 = sin(t365);
t311 = rSges(2,1) * t329 - rSges(2,2) * t326;
t310 = rSges(2,1) * t326 + rSges(2,2) * t329;
t309 = Icges(2,1) * t329 - t374;
t308 = Icges(2,1) * t326 + t323;
t307 = -Icges(2,2) * t326 + t323;
t306 = Icges(2,2) * t329 + t374;
t301 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t300 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t299 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t282 = t377 * rSges(3,3) + (rSges(3,1) * t332 + rSges(3,2) * t380) * t327;
t281 = Icges(3,5) * t377 + (Icges(3,1) * t332 + Icges(3,4) * t380) * t327;
t280 = Icges(3,6) * t377 + (Icges(3,4) * t332 + Icges(3,2) * t380) * t327;
t274 = V_base(5) * rSges(2,3) - t310 * V_base(6) + t364;
t273 = t311 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t269 = t310 * V_base(4) - t311 * V_base(5) + t360;
t268 = -t284 * t352 + t322 * t377;
t267 = -t284 * t322 - t352 * t377;
t255 = rSges(3,1) * t290 + rSges(3,2) * t289 + rSges(3,3) * t373;
t254 = rSges(3,1) * t288 + rSges(3,2) * t287 - rSges(3,3) * t372;
t253 = Icges(3,1) * t290 + Icges(3,4) * t289 + Icges(3,5) * t373;
t252 = Icges(3,1) * t288 + Icges(3,4) * t287 - Icges(3,5) * t372;
t251 = Icges(3,4) * t290 + Icges(3,2) * t289 + Icges(3,6) * t373;
t250 = Icges(3,4) * t288 + Icges(3,2) * t287 - Icges(3,6) * t372;
t247 = -t284 * rSges(4,1) + t283 * rSges(4,2) + rSges(4,3) * t377;
t246 = -Icges(4,1) * t284 + Icges(4,4) * t283 + Icges(4,5) * t377;
t234 = t266 * t352 + t322 * t373;
t233 = t266 * t322 - t326 * t347;
t232 = t264 * t352 - t322 * t372;
t231 = t264 * t322 + t329 * t347;
t230 = t268 * t333 - t283 * t331;
t229 = -t268 * t331 - t283 * t333;
t228 = qJD(6) * t267 + t272;
t227 = pkin(5) * t268 + pkin(9) * t267;
t224 = rSges(5,1) * t271 + rSges(5,2) * t270 - rSges(5,3) * t283;
t223 = Icges(5,1) * t271 + Icges(5,4) * t270 - Icges(5,5) * t283;
t222 = Icges(5,4) * t271 + Icges(5,2) * t270 - Icges(5,6) * t283;
t220 = rSges(4,1) * t266 + rSges(4,2) * t265 + rSges(4,3) * t373;
t219 = rSges(4,1) * t264 + rSges(4,2) * t263 - rSges(4,3) * t372;
t218 = rSges(6,1) * t268 - rSges(6,2) * t267 - rSges(6,3) * t283;
t217 = Icges(4,1) * t266 + Icges(4,4) * t265 + Icges(4,5) * t373;
t216 = Icges(4,1) * t264 + Icges(4,4) * t263 - Icges(4,5) * t372;
t211 = Icges(6,1) * t268 - Icges(6,4) * t267 - Icges(6,5) * t283;
t210 = Icges(6,4) * t268 - Icges(6,2) * t267 - Icges(6,6) * t283;
t209 = Icges(6,5) * t268 - Icges(6,6) * t267 - Icges(6,3) * t283;
t205 = t234 * t333 - t265 * t331;
t204 = -t234 * t331 - t265 * t333;
t203 = t232 * t333 - t263 * t331;
t202 = -t232 * t331 - t263 * t333;
t200 = qJD(6) * t233 + t242;
t199 = qJD(6) * t231 + t241;
t198 = -t254 * t314 + t282 * t302 + t346;
t197 = t314 * t255 - t303 * t282 + t342;
t196 = pkin(5) * t234 + pkin(9) * t233;
t195 = pkin(5) * t232 + pkin(9) * t231;
t194 = t254 * t303 - t255 * t302 + t345;
t193 = rSges(5,1) * t240 + rSges(5,2) * t239 - rSges(5,3) * t265;
t192 = rSges(5,1) * t238 + rSges(5,2) * t237 - rSges(5,3) * t263;
t191 = Icges(5,1) * t240 + Icges(5,4) * t239 - Icges(5,5) * t265;
t190 = Icges(5,1) * t238 + Icges(5,4) * t237 - Icges(5,5) * t263;
t189 = Icges(5,4) * t240 + Icges(5,2) * t239 - Icges(5,6) * t265;
t188 = Icges(5,4) * t238 + Icges(5,2) * t237 - Icges(5,6) * t263;
t185 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t267;
t184 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t267;
t183 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t267;
t182 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t267;
t181 = rSges(6,1) * t234 - rSges(6,2) * t233 - rSges(6,3) * t265;
t180 = rSges(6,1) * t232 - rSges(6,2) * t231 - rSges(6,3) * t263;
t179 = Icges(6,1) * t234 - Icges(6,4) * t233 - Icges(6,5) * t265;
t178 = Icges(6,1) * t232 - Icges(6,4) * t231 - Icges(6,5) * t263;
t177 = Icges(6,4) * t234 - Icges(6,2) * t233 - Icges(6,6) * t265;
t176 = Icges(6,4) * t232 - Icges(6,2) * t231 - Icges(6,6) * t263;
t175 = Icges(6,5) * t234 - Icges(6,6) * t233 - Icges(6,3) * t265;
t174 = Icges(6,5) * t232 - Icges(6,6) * t231 - Icges(6,3) * t263;
t169 = rSges(7,1) * t205 + rSges(7,2) * t204 + rSges(7,3) * t233;
t168 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t231;
t167 = t247 * t302 + (-t219 - t256) * t314 + t344;
t166 = t314 * t220 + (-t247 - t295) * t303 + t338;
t165 = Icges(7,1) * t205 + Icges(7,4) * t204 + Icges(7,5) * t233;
t164 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t231;
t163 = Icges(7,4) * t205 + Icges(7,2) * t204 + Icges(7,6) * t233;
t162 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t231;
t161 = Icges(7,5) * t205 + Icges(7,6) * t204 + Icges(7,3) * t233;
t160 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t231;
t159 = t219 * t303 + (-t220 - t257) * t302 + t343;
t158 = t224 * t302 + (-t192 + t370) * t314 + t341;
t157 = t314 * t193 + (-t224 + t368) * t303 + t337;
t156 = t192 * t303 + (-t193 + t369) * t302 + t340;
t155 = -t180 * t272 + t218 * t241 + t336;
t154 = t272 * t181 - t242 * t218 + t334;
t153 = t180 * t242 - t181 * t241 + t335;
t152 = -t168 * t228 + t185 * t199 - t195 * t272 + t227 * t241 + t336;
t151 = t228 * t169 - t200 * t185 + t272 * t196 - t242 * t227 + t334;
t150 = t168 * t200 - t169 * t199 + t195 * t242 - t196 * t241 + t335;
t1 = m(7) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(6) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(5) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(4) * (t159 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(3) * (t194 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + t199 * ((t161 * t231 + t163 * t202 + t165 * t203) * t200 + (t160 * t231 + t162 * t202 + t164 * t203) * t199 + (t182 * t231 + t183 * t202 + t184 * t203) * t228) / 0.2e1 + t200 * ((t161 * t233 + t163 * t204 + t165 * t205) * t200 + (t160 * t233 + t162 * t204 + t164 * t205) * t199 + (t182 * t233 + t183 * t204 + t184 * t205) * t228) / 0.2e1 + t228 * ((t161 * t267 + t163 * t229 + t165 * t230) * t200 + (t160 * t267 + t162 * t229 + t164 * t230) * t199 + (t182 * t267 + t183 * t229 + t184 * t230) * t228) / 0.2e1 + t241 * ((-t175 * t263 - t177 * t231 + t179 * t232) * t242 + (-t174 * t263 - t176 * t231 + t178 * t232) * t241 + (-t209 * t263 - t210 * t231 + t211 * t232) * t272) / 0.2e1 + t242 * ((-t175 * t265 - t177 * t233 + t179 * t234) * t242 + (-t174 * t265 - t176 * t233 + t178 * t234) * t241 + (-t209 * t265 - t210 * t233 + t211 * t234) * t272) / 0.2e1 + m(2) * (t269 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + t272 * ((-t175 * t283 - t177 * t267 + t179 * t268) * t242 + (-t174 * t283 - t176 * t267 + t178 * t268) * t241 + (-t209 * t283 - t210 * t267 + t211 * t268) * t272) / 0.2e1 + m(1) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + ((-t306 * t326 + t308 * t329 + Icges(1,4)) * V_base(5) + (-t307 * t326 + t309 * t329 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t306 * t329 + t308 * t326 + Icges(1,2)) * V_base(5) + (t307 * t329 + t309 * t326 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t222 * t237 + t223 * t238 + t246 * t264 - t263 * t384 + t280 * t287 + t281 * t288 - t372 * t385) * t314 + (t189 * t237 + t191 * t238 + t217 * t264 + t251 * t287 + t253 * t288 - t263 * t386 - t372 * t388) * t303 + (t188 * t237 + t190 * t238 + t216 * t264 + t250 * t287 + t252 * t288 - t387 * t263 - t389 * t372) * t302) * t302 / 0.2e1 + ((t222 * t239 + t223 * t240 + t246 * t266 - t265 * t384 + t280 * t289 + t281 * t290 + t373 * t385) * t314 + (t189 * t239 + t191 * t240 + t217 * t266 + t251 * t289 + t253 * t290 - t386 * t265 + t388 * t373) * t303 + (t188 * t239 + t190 * t240 + t216 * t266 + t250 * t289 + t252 * t290 - t265 * t387 + t373 * t389) * t302) * t303 / 0.2e1 + ((-t284 * t246 + t222 * t270 + t223 * t271 + (t280 * t380 + t281 * t332) * t327 + t385 * t377 - t384 * t283) * t314 + (-t284 * t217 + t189 * t270 + t191 * t271 + (t251 * t380 + t253 * t332) * t327 + t388 * t377 - t386 * t283) * t303 + (-t284 * t216 + t188 * t270 + t190 * t271 + (t250 * t380 + t252 * t332) * t327 + t389 * t377 - t387 * t283) * t302) * t314 / 0.2e1 + ((Icges(2,5) * t326 + Icges(2,6) * t329 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t329 - Icges(2,6) * t326 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
