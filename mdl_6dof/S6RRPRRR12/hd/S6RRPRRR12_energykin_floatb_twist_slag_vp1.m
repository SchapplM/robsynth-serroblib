% Calculate kinetic energy for
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:08
% EndTime: 2019-03-09 14:35:12
% DurationCPUTime: 4.44s
% Computational Cost: add. (2647->396), mult. (4786->574), div. (0->0), fcn. (5502->12), ass. (0->177)
t366 = Icges(3,1) + Icges(4,2);
t365 = Icges(3,4) + Icges(4,6);
t364 = Icges(3,5) - Icges(4,4);
t363 = Icges(3,2) + Icges(4,3);
t362 = Icges(3,6) - Icges(4,5);
t361 = Icges(3,3) + Icges(4,1);
t305 = cos(pkin(6));
t309 = sin(qJ(1));
t312 = cos(qJ(2));
t335 = t309 * t312;
t308 = sin(qJ(2));
t313 = cos(qJ(1));
t336 = t308 * t313;
t271 = t305 * t335 + t336;
t334 = t312 * t313;
t337 = t308 * t309;
t272 = -t305 * t337 + t334;
t304 = sin(pkin(6));
t341 = t304 * t309;
t360 = t363 * t271 - t365 * t272 - t362 * t341;
t269 = -t305 * t334 + t337;
t270 = t305 * t336 + t335;
t339 = t304 * t313;
t359 = t363 * t269 - t365 * t270 + t362 * t339;
t358 = -t365 * t271 + t366 * t272 + t364 * t341;
t357 = -t365 * t269 + t366 * t270 - t364 * t339;
t356 = -t362 * t271 + t364 * t272 + t361 * t341;
t355 = -t362 * t269 + t364 * t270 - t361 * t339;
t354 = t361 * t305 + (t364 * t308 + t362 * t312) * t304;
t353 = t362 * t305 + (t365 * t308 + t363 * t312) * t304;
t352 = t364 * t305 + (t366 * t308 + t365 * t312) * t304;
t348 = pkin(8) * t305;
t311 = cos(qJ(4));
t347 = pkin(4) * t311;
t345 = Icges(2,4) * t309;
t307 = sin(qJ(4));
t344 = t269 * t307;
t343 = t271 * t307;
t342 = t304 * t308;
t340 = t304 * t312;
t338 = t307 * t312;
t333 = qJ(4) + qJ(5);
t332 = qJD(2) * t304;
t331 = V_base(5) * pkin(7) + V_base(1);
t282 = t309 * t332 + V_base(4);
t328 = cos(t333);
t300 = V_base(6) + qJD(1);
t235 = qJD(4) * t272 + t282;
t283 = qJD(2) * t305 + t300;
t327 = t304 * t328;
t200 = qJD(5) * t272 + t235;
t265 = qJD(4) * t342 + t283;
t281 = -t313 * t332 + V_base(5);
t276 = pkin(1) * t309 - pkin(8) * t339;
t326 = -t276 * t300 + V_base(5) * t348 + t331;
t244 = qJD(5) * t342 + t265;
t277 = pkin(1) * t313 + pkin(8) * t341;
t325 = V_base(4) * t276 - t277 * V_base(5) + V_base(3);
t234 = qJD(4) * t270 + t281;
t199 = qJD(5) * t270 + t234;
t273 = (pkin(2) * t308 - qJ(3) * t312) * t304;
t324 = qJD(3) * t271 + t281 * t273 + t326;
t323 = t300 * t277 + V_base(2) + (-pkin(7) - t348) * V_base(4);
t227 = pkin(2) * t272 + qJ(3) * t271;
t322 = qJD(3) * t269 + t283 * t227 + t323;
t226 = pkin(2) * t270 + qJ(3) * t269;
t321 = -qJD(3) * t340 + t282 * t226 + t325;
t246 = -pkin(3) * t339 + t270 * pkin(9);
t275 = t305 * pkin(3) + pkin(9) * t342;
t320 = t281 * t275 + (-t226 - t246) * t283 + t324;
t245 = pkin(3) * t341 + t272 * pkin(9);
t319 = t283 * t245 + (-t273 - t275) * t282 + t322;
t318 = t282 * t246 + (-t227 - t245) * t281 + t321;
t185 = pkin(4) * t344 + pkin(10) * t270 - t339 * t347;
t223 = t347 * t305 + (-pkin(4) * t338 + pkin(10) * t308) * t304;
t317 = -t185 * t265 + t234 * t223 + t320;
t184 = pkin(4) * t343 + pkin(10) * t272 + t341 * t347;
t316 = t265 * t184 - t223 * t235 + t319;
t315 = -t184 * t234 + t235 * t185 + t318;
t310 = cos(qJ(6));
t306 = sin(qJ(6));
t302 = Icges(2,4) * t313;
t301 = sin(t333);
t291 = rSges(2,1) * t313 - rSges(2,2) * t309;
t290 = rSges(2,1) * t309 + rSges(2,2) * t313;
t289 = Icges(2,1) * t313 - t345;
t288 = Icges(2,1) * t309 + t302;
t287 = -Icges(2,2) * t309 + t302;
t286 = Icges(2,2) * t313 + t345;
t280 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t279 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t278 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t268 = -t304 * t338 + t305 * t311;
t267 = -t305 * t307 - t311 * t340;
t258 = -t301 * t340 + t305 * t328;
t257 = t305 * t301 + t312 * t327;
t256 = rSges(4,1) * t305 + (-rSges(4,2) * t308 - rSges(4,3) * t312) * t304;
t255 = rSges(3,3) * t305 + (rSges(3,1) * t308 + rSges(3,2) * t312) * t304;
t243 = V_base(5) * rSges(2,3) - t290 * t300 + t331;
t242 = t291 * t300 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t240 = t290 * V_base(4) - t291 * V_base(5) + V_base(3);
t239 = -t311 * t339 + t344;
t238 = t269 * t311 + t307 * t339;
t237 = t311 * t341 + t343;
t236 = t271 * t311 - t307 * t341;
t233 = t269 * t301 - t313 * t327;
t232 = t269 * t328 + t301 * t339;
t231 = t271 * t301 + t309 * t327;
t230 = -t271 * t328 + t301 * t341;
t229 = t258 * t310 + t306 * t342;
t228 = -t258 * t306 + t310 * t342;
t222 = rSges(3,1) * t272 - rSges(3,2) * t271 + rSges(3,3) * t341;
t221 = rSges(3,1) * t270 - rSges(3,2) * t269 - rSges(3,3) * t339;
t220 = -rSges(4,1) * t339 - rSges(4,2) * t270 + rSges(4,3) * t269;
t219 = rSges(4,1) * t341 - rSges(4,2) * t272 + rSges(4,3) * t271;
t206 = pkin(5) * t258 + pkin(11) * t257;
t205 = rSges(5,1) * t268 + rSges(5,2) * t267 + rSges(5,3) * t342;
t204 = Icges(5,1) * t268 + Icges(5,4) * t267 + Icges(5,5) * t342;
t203 = Icges(5,4) * t268 + Icges(5,2) * t267 + Icges(5,6) * t342;
t202 = Icges(5,5) * t268 + Icges(5,6) * t267 + Icges(5,3) * t342;
t197 = qJD(6) * t257 + t244;
t196 = rSges(6,1) * t258 - rSges(6,2) * t257 + rSges(6,3) * t342;
t195 = Icges(6,1) * t258 - Icges(6,4) * t257 + Icges(6,5) * t342;
t194 = Icges(6,4) * t258 - Icges(6,2) * t257 + Icges(6,6) * t342;
t193 = Icges(6,5) * t258 - Icges(6,6) * t257 + Icges(6,3) * t342;
t192 = t233 * t310 + t270 * t306;
t191 = -t233 * t306 + t270 * t310;
t190 = t231 * t310 + t272 * t306;
t189 = -t231 * t306 + t272 * t310;
t188 = pkin(5) * t233 - pkin(11) * t232;
t187 = pkin(5) * t231 + pkin(11) * t230;
t183 = rSges(5,1) * t239 + rSges(5,2) * t238 + rSges(5,3) * t270;
t182 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t272;
t181 = Icges(5,1) * t239 + Icges(5,4) * t238 + Icges(5,5) * t270;
t180 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t272;
t179 = Icges(5,4) * t239 + Icges(5,2) * t238 + Icges(5,6) * t270;
t178 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t272;
t177 = Icges(5,5) * t239 + Icges(5,6) * t238 + Icges(5,3) * t270;
t176 = Icges(5,5) * t237 + Icges(5,6) * t236 + Icges(5,3) * t272;
t175 = qJD(6) * t230 + t200;
t174 = -qJD(6) * t232 + t199;
t173 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t270;
t172 = rSges(6,1) * t231 - rSges(6,2) * t230 + rSges(6,3) * t272;
t171 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t270;
t170 = Icges(6,1) * t231 - Icges(6,4) * t230 + Icges(6,5) * t272;
t169 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t270;
t168 = Icges(6,4) * t231 - Icges(6,2) * t230 + Icges(6,6) * t272;
t167 = Icges(6,5) * t233 + Icges(6,6) * t232 + Icges(6,3) * t270;
t166 = Icges(6,5) * t231 - Icges(6,6) * t230 + Icges(6,3) * t272;
t165 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t257;
t164 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t257;
t163 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t257;
t162 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t257;
t159 = -t221 * t283 + t255 * t281 + t326;
t158 = t222 * t283 - t255 * t282 + t323;
t157 = t221 * t282 - t222 * t281 + t325;
t156 = rSges(7,1) * t192 + rSges(7,2) * t191 - rSges(7,3) * t232;
t155 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t230;
t154 = Icges(7,1) * t192 + Icges(7,4) * t191 - Icges(7,5) * t232;
t153 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t230;
t152 = Icges(7,4) * t192 + Icges(7,2) * t191 - Icges(7,6) * t232;
t151 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t230;
t150 = Icges(7,5) * t192 + Icges(7,6) * t191 - Icges(7,3) * t232;
t149 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t230;
t148 = t256 * t281 + (-t220 - t226) * t283 + t324;
t147 = t219 * t283 + (-t256 - t273) * t282 + t322;
t146 = t220 * t282 + (-t219 - t227) * t281 + t321;
t145 = -t183 * t265 + t205 * t234 + t320;
t144 = t182 * t265 - t205 * t235 + t319;
t143 = -t182 * t234 + t183 * t235 + t318;
t142 = -t173 * t244 + t196 * t199 + t317;
t141 = t172 * t244 - t196 * t200 + t316;
t140 = -t172 * t199 + t173 * t200 + t315;
t139 = -t156 * t197 + t165 * t174 - t188 * t244 + t199 * t206 + t317;
t138 = t155 * t197 - t165 * t175 + t187 * t244 - t200 * t206 + t316;
t137 = -t155 * t174 + t156 * t175 - t187 * t199 + t188 * t200 + t315;
t1 = t265 * ((t176 * t342 + t178 * t267 + t180 * t268) * t235 + (t177 * t342 + t179 * t267 + t181 * t268) * t234 + (t202 * t342 + t267 * t203 + t268 * t204) * t265) / 0.2e1 + t244 * ((t166 * t342 - t168 * t257 + t170 * t258) * t200 + (t167 * t342 - t169 * t257 + t171 * t258) * t199 + (t193 * t342 - t257 * t194 + t258 * t195) * t244) / 0.2e1 + m(1) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + t235 * ((t272 * t176 + t236 * t178 + t237 * t180) * t235 + (t177 * t272 + t179 * t236 + t181 * t237) * t234 + (t202 * t272 + t203 * t236 + t204 * t237) * t265) / 0.2e1 + t200 * ((t272 * t166 - t230 * t168 + t231 * t170) * t200 + (t167 * t272 - t169 * t230 + t171 * t231) * t199 + (t193 * t272 - t194 * t230 + t195 * t231) * t244) / 0.2e1 + t234 * ((t176 * t270 + t178 * t238 + t180 * t239) * t235 + (t270 * t177 + t238 * t179 + t239 * t181) * t234 + (t202 * t270 + t203 * t238 + t204 * t239) * t265) / 0.2e1 + t199 * ((t166 * t270 + t168 * t232 + t170 * t233) * t200 + (t270 * t167 + t232 * t169 + t233 * t171) * t199 + (t193 * t270 + t194 * t232 + t195 * t233) * t244) / 0.2e1 + t197 * ((t149 * t257 + t151 * t228 + t153 * t229) * t175 + (t150 * t257 + t152 * t228 + t154 * t229) * t174 + (t257 * t162 + t228 * t163 + t229 * t164) * t197) / 0.2e1 + m(2) * (t240 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + t175 * ((t230 * t149 + t189 * t151 + t190 * t153) * t175 + (t150 * t230 + t152 * t189 + t154 * t190) * t174 + (t162 * t230 + t163 * t189 + t164 * t190) * t197) / 0.2e1 + t174 * ((-t149 * t232 + t151 * t191 + t153 * t192) * t175 + (-t232 * t150 + t191 * t152 + t192 * t154) * t174 + (-t162 * t232 + t163 * t191 + t164 * t192) * t197) / 0.2e1 + m(3) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(5) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(4) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(6) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(7) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + ((-t269 * t353 + t270 * t352 - t339 * t354) * t283 + (t269 * t360 + t358 * t270 - t356 * t339) * t282 + (t359 * t269 + t357 * t270 - t355 * t339) * t281) * t281 / 0.2e1 + ((-t271 * t353 + t272 * t352 + t341 * t354) * t283 + (t360 * t271 + t358 * t272 + t356 * t341) * t282 + (t271 * t359 + t272 * t357 + t341 * t355) * t281) * t282 / 0.2e1 + ((t281 * t355 + t282 * t356 + t283 * t354) * t305 + ((t308 * t352 + t312 * t353) * t283 + (t358 * t308 - t312 * t360) * t282 + (t308 * t357 - t312 * t359) * t281) * t304) * t283 / 0.2e1 + ((-t286 * t309 + t288 * t313 + Icges(1,4)) * V_base(5) + (-t287 * t309 + t289 * t313 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t286 * t313 + t288 * t309 + Icges(1,2)) * V_base(5) + (t287 * t313 + t289 * t309 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t309 + Icges(2,6) * t313) * V_base(5) + (Icges(2,5) * t313 - Icges(2,6) * t309) * V_base(4) + Icges(2,3) * t300 / 0.2e1) * t300;
T  = t1;
