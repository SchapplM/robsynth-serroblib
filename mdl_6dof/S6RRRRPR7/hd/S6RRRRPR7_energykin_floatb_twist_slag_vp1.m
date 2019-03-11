% Calculate kinetic energy for
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:26
% EndTime: 2019-03-09 22:26:30
% DurationCPUTime: 4.80s
% Computational Cost: add. (3737->438), mult. (5584->633), div. (0->0), fcn. (6496->14), ass. (0->190)
t373 = Icges(5,3) + Icges(6,3);
t318 = cos(pkin(6));
t322 = sin(qJ(1));
t325 = cos(qJ(2));
t353 = t322 * t325;
t321 = sin(qJ(2));
t326 = cos(qJ(1));
t354 = t321 * t326;
t278 = t318 * t354 + t353;
t316 = qJ(3) + qJ(4);
t343 = pkin(12) + t316;
t306 = sin(t343);
t317 = sin(pkin(6));
t341 = cos(t343);
t340 = t317 * t341;
t238 = t278 * t306 + t326 * t340;
t357 = t317 * t326;
t239 = t278 * t341 - t306 * t357;
t310 = sin(t316);
t311 = cos(t316);
t244 = -t278 * t310 - t311 * t357;
t245 = t278 * t311 - t310 * t357;
t352 = t325 * t326;
t355 = t321 * t322;
t277 = -t318 * t352 + t355;
t372 = Icges(5,5) * t245 + Icges(6,5) * t239 + Icges(5,6) * t244 - Icges(6,6) * t238 + t373 * t277;
t280 = -t318 * t355 + t352;
t240 = t280 * t306 - t322 * t340;
t360 = t317 * t322;
t241 = t280 * t341 + t306 * t360;
t246 = -t280 * t310 + t311 * t360;
t247 = t280 * t311 + t310 * t360;
t279 = t318 * t353 + t354;
t371 = Icges(5,5) * t247 + Icges(6,5) * t241 + Icges(5,6) * t246 - Icges(6,6) * t240 + t373 * t279;
t361 = t317 * t321;
t259 = t306 * t361 - t318 * t341;
t260 = t318 * t306 + t321 * t340;
t265 = -t310 * t361 + t311 * t318;
t266 = t310 * t318 + t311 * t361;
t358 = t317 * t325;
t370 = Icges(5,5) * t266 + Icges(6,5) * t260 + Icges(5,6) * t265 - Icges(6,6) * t259 - t373 * t358;
t365 = pkin(8) * t318;
t324 = cos(qJ(3));
t364 = t324 * pkin(3);
t362 = Icges(2,4) * t322;
t359 = t317 * t324;
t320 = sin(qJ(3));
t356 = t318 * t320;
t351 = pkin(4) * t311;
t349 = qJD(2) * t317;
t348 = V_base(5) * pkin(7) + V_base(1);
t345 = t320 * t360;
t344 = t320 * t357;
t292 = t322 * t349 + V_base(4);
t309 = V_base(6) + qJD(1);
t342 = pkin(4) * t310;
t249 = qJD(3) * t279 + t292;
t293 = qJD(2) * t318 + t309;
t221 = qJD(4) * t279 + t249;
t291 = -t326 * t349 + V_base(5);
t284 = pkin(1) * t322 - pkin(8) * t357;
t339 = -t284 * t309 + V_base(5) * t365 + t348;
t285 = pkin(1) * t326 + pkin(8) * t360;
t338 = V_base(4) * t284 - t285 * V_base(5) + V_base(3);
t248 = qJD(3) * t277 + t291;
t220 = qJD(4) * t277 + t248;
t337 = t309 * t285 + V_base(2) + (-pkin(7) - t365) * V_base(4);
t258 = (-qJD(3) - qJD(4)) * t358 + t293;
t242 = t278 * pkin(2) + t277 * pkin(9);
t282 = (pkin(2) * t321 - pkin(9) * t325) * t317;
t336 = -t242 * t293 + t291 * t282 + t339;
t243 = t280 * pkin(2) + t279 * pkin(9);
t335 = t292 * t242 - t243 * t291 + t338;
t334 = t293 * t243 - t282 * t292 + t337;
t189 = -pkin(3) * t344 + pkin(10) * t277 + t278 * t364;
t233 = pkin(3) * t356 + (-pkin(10) * t325 + t321 * t364) * t317;
t273 = -qJD(3) * t358 + t293;
t333 = -t189 * t273 + t248 * t233 + t336;
t190 = pkin(3) * t345 + pkin(10) * t279 + t280 * t364;
t332 = t249 * t189 - t190 * t248 + t335;
t204 = t342 * t318 + (-qJ(5) * t325 + t321 * t351) * t317;
t331 = qJD(5) * t279 + t220 * t204 + t333;
t330 = t273 * t190 - t233 * t249 + t334;
t164 = qJ(5) * t279 + t280 * t351 + t342 * t360;
t329 = qJD(5) * t277 + t258 * t164 + t330;
t163 = qJ(5) * t277 + t278 * t351 - t342 * t357;
t328 = -qJD(5) * t358 + t221 * t163 + t332;
t323 = cos(qJ(6));
t319 = sin(qJ(6));
t312 = Icges(2,4) * t326;
t302 = rSges(2,1) * t326 - rSges(2,2) * t322;
t301 = rSges(2,1) * t322 + rSges(2,2) * t326;
t299 = Icges(2,1) * t326 - t362;
t298 = Icges(2,1) * t322 + t312;
t297 = -Icges(2,2) * t322 + t312;
t296 = Icges(2,2) * t326 + t362;
t289 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t288 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t287 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t276 = t321 * t359 + t356;
t275 = t318 * t324 - t320 * t361;
t264 = rSges(3,3) * t318 + (rSges(3,1) * t321 + rSges(3,2) * t325) * t317;
t263 = Icges(3,5) * t318 + (Icges(3,1) * t321 + Icges(3,4) * t325) * t317;
t262 = Icges(3,6) * t318 + (Icges(3,4) * t321 + Icges(3,2) * t325) * t317;
t261 = Icges(3,3) * t318 + (Icges(3,5) * t321 + Icges(3,6) * t325) * t317;
t257 = V_base(5) * rSges(2,3) - t301 * t309 + t348;
t256 = t302 * t309 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t254 = t301 * V_base(4) - t302 * V_base(5) + V_base(3);
t253 = t280 * t324 + t345;
t252 = -t280 * t320 + t322 * t359;
t251 = t278 * t324 - t344;
t250 = -t278 * t320 - t324 * t357;
t237 = t260 * t323 - t319 * t358;
t236 = -t260 * t319 - t323 * t358;
t235 = rSges(3,1) * t280 - rSges(3,2) * t279 + rSges(3,3) * t360;
t234 = rSges(3,1) * t278 - rSges(3,2) * t277 - rSges(3,3) * t357;
t232 = Icges(3,1) * t280 - Icges(3,4) * t279 + Icges(3,5) * t360;
t231 = Icges(3,1) * t278 - Icges(3,4) * t277 - Icges(3,5) * t357;
t230 = Icges(3,4) * t280 - Icges(3,2) * t279 + Icges(3,6) * t360;
t229 = Icges(3,4) * t278 - Icges(3,2) * t277 - Icges(3,6) * t357;
t228 = Icges(3,5) * t280 - Icges(3,6) * t279 + Icges(3,3) * t360;
t227 = Icges(3,5) * t278 - Icges(3,6) * t277 - Icges(3,3) * t357;
t226 = rSges(4,1) * t276 + rSges(4,2) * t275 - rSges(4,3) * t358;
t225 = Icges(4,1) * t276 + Icges(4,4) * t275 - Icges(4,5) * t358;
t224 = Icges(4,4) * t276 + Icges(4,2) * t275 - Icges(4,6) * t358;
t223 = Icges(4,5) * t276 + Icges(4,6) * t275 - Icges(4,3) * t358;
t218 = pkin(5) * t260 + pkin(11) * t259;
t217 = qJD(6) * t259 + t258;
t216 = rSges(5,1) * t266 + rSges(5,2) * t265 - rSges(5,3) * t358;
t215 = Icges(5,1) * t266 + Icges(5,4) * t265 - Icges(5,5) * t358;
t214 = Icges(5,4) * t266 + Icges(5,2) * t265 - Icges(5,6) * t358;
t212 = rSges(6,1) * t260 - rSges(6,2) * t259 - rSges(6,3) * t358;
t211 = Icges(6,1) * t260 - Icges(6,4) * t259 - Icges(6,5) * t358;
t210 = Icges(6,4) * t260 - Icges(6,2) * t259 - Icges(6,6) * t358;
t208 = t241 * t323 + t279 * t319;
t207 = -t241 * t319 + t279 * t323;
t206 = t239 * t323 + t277 * t319;
t205 = -t239 * t319 + t277 * t323;
t202 = pkin(5) * t241 + pkin(11) * t240;
t201 = pkin(5) * t239 + pkin(11) * t238;
t200 = rSges(4,1) * t253 + rSges(4,2) * t252 + rSges(4,3) * t279;
t199 = rSges(4,1) * t251 + rSges(4,2) * t250 + rSges(4,3) * t277;
t198 = Icges(4,1) * t253 + Icges(4,4) * t252 + Icges(4,5) * t279;
t197 = Icges(4,1) * t251 + Icges(4,4) * t250 + Icges(4,5) * t277;
t196 = Icges(4,4) * t253 + Icges(4,2) * t252 + Icges(4,6) * t279;
t195 = Icges(4,4) * t251 + Icges(4,2) * t250 + Icges(4,6) * t277;
t194 = Icges(4,5) * t253 + Icges(4,6) * t252 + Icges(4,3) * t279;
t193 = Icges(4,5) * t251 + Icges(4,6) * t250 + Icges(4,3) * t277;
t192 = qJD(6) * t240 + t221;
t191 = qJD(6) * t238 + t220;
t188 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t279;
t187 = rSges(5,1) * t245 + rSges(5,2) * t244 + rSges(5,3) * t277;
t186 = Icges(5,1) * t247 + Icges(5,4) * t246 + Icges(5,5) * t279;
t185 = Icges(5,1) * t245 + Icges(5,4) * t244 + Icges(5,5) * t277;
t184 = Icges(5,4) * t247 + Icges(5,2) * t246 + Icges(5,6) * t279;
t183 = Icges(5,4) * t245 + Icges(5,2) * t244 + Icges(5,6) * t277;
t180 = rSges(6,1) * t241 - rSges(6,2) * t240 + rSges(6,3) * t279;
t179 = rSges(6,1) * t239 - rSges(6,2) * t238 + rSges(6,3) * t277;
t178 = Icges(6,1) * t241 - Icges(6,4) * t240 + Icges(6,5) * t279;
t177 = Icges(6,1) * t239 - Icges(6,4) * t238 + Icges(6,5) * t277;
t176 = Icges(6,4) * t241 - Icges(6,2) * t240 + Icges(6,6) * t279;
t175 = Icges(6,4) * t239 - Icges(6,2) * t238 + Icges(6,6) * t277;
t172 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t259;
t171 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t259;
t170 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t259;
t169 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t259;
t166 = -t234 * t293 + t264 * t291 + t339;
t165 = t235 * t293 - t264 * t292 + t337;
t161 = t234 * t292 - t235 * t291 + t338;
t159 = rSges(7,1) * t208 + rSges(7,2) * t207 + rSges(7,3) * t240;
t158 = rSges(7,1) * t206 + rSges(7,2) * t205 + rSges(7,3) * t238;
t157 = Icges(7,1) * t208 + Icges(7,4) * t207 + Icges(7,5) * t240;
t156 = Icges(7,1) * t206 + Icges(7,4) * t205 + Icges(7,5) * t238;
t155 = Icges(7,4) * t208 + Icges(7,2) * t207 + Icges(7,6) * t240;
t154 = Icges(7,4) * t206 + Icges(7,2) * t205 + Icges(7,6) * t238;
t153 = Icges(7,5) * t208 + Icges(7,6) * t207 + Icges(7,3) * t240;
t152 = Icges(7,5) * t206 + Icges(7,6) * t205 + Icges(7,3) * t238;
t150 = -t199 * t273 + t226 * t248 + t336;
t149 = t200 * t273 - t226 * t249 + t334;
t148 = t199 * t249 - t200 * t248 + t335;
t147 = -t187 * t258 + t216 * t220 + t333;
t146 = t188 * t258 - t216 * t221 + t330;
t145 = t187 * t221 - t188 * t220 + t332;
t144 = t212 * t220 + (-t163 - t179) * t258 + t331;
t143 = t180 * t258 + (-t204 - t212) * t221 + t329;
t142 = t179 * t221 + (-t164 - t180) * t220 + t328;
t141 = (-t163 - t201) * t258 - t158 * t217 + t172 * t191 + t218 * t220 + t331;
t140 = t159 * t217 - t172 * t192 + t202 * t258 + (-t204 - t218) * t221 + t329;
t139 = t158 * t192 - t159 * t191 + t201 * t221 + (-t164 - t202) * t220 + t328;
t1 = t292 * ((t228 * t360 - t230 * t279 + t232 * t280) * t292 + (t227 * t360 - t229 * t279 + t231 * t280) * t291 + (t261 * t360 - t262 * t279 + t263 * t280) * t293) / 0.2e1 + t273 * ((-t194 * t358 + t196 * t275 + t198 * t276) * t249 + (-t193 * t358 + t195 * t275 + t197 * t276) * t248 + (-t223 * t358 + t224 * t275 + t225 * t276) * t273) / 0.2e1 + t291 * ((-t228 * t357 - t230 * t277 + t232 * t278) * t292 + (-t227 * t357 - t229 * t277 + t231 * t278) * t291 + (-t261 * t357 - t262 * t277 + t263 * t278) * t293) / 0.2e1 + m(1) * (t287 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + t249 * ((t194 * t279 + t196 * t252 + t198 * t253) * t249 + (t193 * t279 + t195 * t252 + t197 * t253) * t248 + (t223 * t279 + t224 * t252 + t225 * t253) * t273) / 0.2e1 + t248 * ((t194 * t277 + t196 * t250 + t198 * t251) * t249 + (t193 * t277 + t195 * t250 + t197 * t251) * t248 + (t223 * t277 + t224 * t250 + t225 * t251) * t273) / 0.2e1 + m(2) * (t254 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t217 * ((t153 * t259 + t155 * t236 + t157 * t237) * t192 + (t152 * t259 + t154 * t236 + t156 * t237) * t191 + (t169 * t259 + t170 * t236 + t171 * t237) * t217) / 0.2e1 + t191 * ((t153 * t238 + t155 * t205 + t157 * t206) * t192 + (t152 * t238 + t154 * t205 + t156 * t206) * t191 + (t169 * t238 + t170 * t205 + t171 * t206) * t217) / 0.2e1 + t192 * ((t153 * t240 + t155 * t207 + t157 * t208) * t192 + (t152 * t240 + t154 * t207 + t156 * t208) * t191 + (t169 * t240 + t170 * t207 + t171 * t208) * t217) / 0.2e1 + m(3) * (t161 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + t293 * ((t227 * t291 + t228 * t292 + t261 * t293) * t318 + ((t230 * t325 + t232 * t321) * t292 + (t229 * t325 + t231 * t321) * t291 + (t262 * t325 + t263 * t321) * t293) * t317) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + ((-t210 * t238 + t211 * t239 + t214 * t244 + t215 * t245 + t277 * t370) * t258 + (-t176 * t238 + t178 * t239 + t184 * t244 + t186 * t245 + t277 * t371) * t221 + (-t175 * t238 + t177 * t239 + t183 * t244 + t185 * t245 + t372 * t277) * t220) * t220 / 0.2e1 + ((-t210 * t240 + t211 * t241 + t214 * t246 + t215 * t247 + t279 * t370) * t258 + (-t176 * t240 + t178 * t241 + t184 * t246 + t186 * t247 + t371 * t279) * t221 + (-t175 * t240 + t177 * t241 + t183 * t246 + t185 * t247 + t279 * t372) * t220) * t221 / 0.2e1 + ((-t210 * t259 + t211 * t260 + t214 * t265 + t215 * t266 - t370 * t358) * t258 + (-t176 * t259 + t178 * t260 + t184 * t265 + t186 * t266 - t358 * t371) * t221 + (-t175 * t259 + t177 * t260 + t183 * t265 + t185 * t266 - t358 * t372) * t220) * t258 / 0.2e1 + ((-t296 * t322 + t298 * t326 + Icges(1,4)) * V_base(5) + (-t297 * t322 + t299 * t326 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t326 + t298 * t322 + Icges(1,2)) * V_base(5) + (t297 * t326 + t299 * t322 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t322 + Icges(2,6) * t326) * V_base(5) + (Icges(2,5) * t326 - Icges(2,6) * t322) * V_base(4) + Icges(2,3) * t309 / 0.2e1) * t309;
T  = t1;
