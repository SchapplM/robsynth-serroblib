% Calculate kinetic energy for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:14
% EndTime: 2019-03-09 11:47:18
% DurationCPUTime: 3.50s
% Computational Cost: add. (2301->335), mult. (2353->480), div. (0->0), fcn. (2247->10), ass. (0->168)
t349 = Icges(6,1) + Icges(7,1);
t348 = Icges(6,4) + Icges(7,4);
t347 = -Icges(7,5) - Icges(6,5);
t346 = Icges(6,2) + Icges(7,2);
t345 = -Icges(7,6) - Icges(6,6);
t344 = Icges(3,3) + Icges(4,3);
t343 = -Icges(7,3) - Icges(6,3);
t255 = qJ(2) + pkin(10);
t244 = sin(t255);
t245 = cos(t255);
t259 = sin(qJ(2));
t262 = cos(qJ(2));
t342 = Icges(3,5) * t262 + Icges(4,5) * t245 - Icges(3,6) * t259 - Icges(4,6) * t244;
t256 = qJ(4) + qJ(5);
t250 = cos(t256);
t263 = cos(qJ(1));
t308 = t250 * t263;
t249 = sin(t256);
t260 = sin(qJ(1));
t311 = t249 * t260;
t189 = -t245 * t311 - t308;
t309 = t250 * t260;
t310 = t249 * t263;
t190 = t245 * t309 - t310;
t313 = t244 * t260;
t341 = -t345 * t189 - t347 * t190 - t343 * t313;
t191 = -t245 * t310 + t309;
t192 = t245 * t308 + t311;
t312 = t244 * t263;
t340 = -t345 * t191 - t347 * t192 - t343 * t312;
t339 = t346 * t189 + t348 * t190 - t345 * t313;
t338 = t346 * t191 + t348 * t192 - t345 * t312;
t337 = t348 * t189 + t349 * t190 - t347 * t313;
t336 = t348 * t191 + t349 * t192 - t347 * t312;
t335 = t343 * t245 + (t345 * t249 - t347 * t250) * t244;
t334 = t345 * t245 + (-t346 * t249 + t348 * t250) * t244;
t333 = t347 * t245 + (-t348 * t249 + t349 * t250) * t244;
t314 = Icges(4,4) * t245;
t283 = -Icges(4,2) * t244 + t314;
t182 = -Icges(4,6) * t263 + t260 * t283;
t183 = Icges(4,6) * t260 + t263 * t283;
t315 = Icges(4,4) * t244;
t285 = Icges(4,1) * t245 - t315;
t184 = -Icges(4,5) * t263 + t260 * t285;
t185 = Icges(4,5) * t260 + t263 * t285;
t316 = Icges(3,4) * t262;
t284 = -Icges(3,2) * t259 + t316;
t195 = -Icges(3,6) * t263 + t260 * t284;
t196 = Icges(3,6) * t260 + t263 * t284;
t317 = Icges(3,4) * t259;
t286 = Icges(3,1) * t262 - t317;
t197 = -Icges(3,5) * t263 + t260 * t286;
t198 = Icges(3,5) * t260 + t263 * t286;
t212 = Icges(4,2) * t245 + t315;
t213 = Icges(4,1) * t244 + t314;
t227 = Icges(3,2) * t262 + t317;
t230 = Icges(3,1) * t259 + t316;
t240 = -qJD(2) * t263 + V_base(5);
t241 = qJD(2) * t260 + V_base(4);
t246 = V_base(6) + qJD(1);
t332 = (-t212 * t244 + t213 * t245 - t227 * t259 + t230 * t262) * t246 + (-t183 * t244 + t185 * t245 - t196 * t259 + t198 * t262) * t241 + (-t182 * t244 + t184 * t245 - t195 * t259 + t197 * t262) * t240;
t331 = (Icges(3,5) * t259 + Icges(4,5) * t244 + Icges(3,6) * t262 + Icges(4,6) * t245) * t246 + (t344 * t260 + t342 * t263) * t241 + (t342 * t260 - t344 * t263) * t240;
t324 = pkin(2) * t259;
t322 = pkin(2) * t262;
t261 = cos(qJ(4));
t321 = t261 * pkin(4);
t318 = Icges(2,4) * t260;
t258 = sin(qJ(4));
t307 = t258 * t260;
t306 = t258 * t263;
t305 = t260 * t261;
t304 = t261 * t263;
t299 = pkin(5) * t250;
t274 = qJ(6) * t244 + t245 * t299;
t291 = pkin(5) * t249;
t303 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t313 + t260 * t274 - t263 * t291;
t302 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t312 + t260 * t291 + t263 * t274;
t301 = (-qJ(6) - rSges(7,3)) * t245 + (rSges(7,1) * t250 - rSges(7,2) * t249 + t299) * t244;
t177 = -qJ(3) * t263 + t260 * t322;
t238 = pkin(1) * t260 - pkin(7) * t263;
t300 = -t177 - t238;
t297 = qJD(4) * t244;
t296 = qJD(5) * t244;
t295 = qJD(6) * t244;
t294 = V_base(5) * pkin(6) + V_base(1);
t205 = t263 * t297 + t241;
t290 = qJD(3) * t260 + t240 * t324 + t294;
t289 = pkin(3) * t245 + pkin(8) * t244;
t288 = rSges(3,1) * t262 - rSges(3,2) * t259;
t287 = rSges(4,1) * t245 - rSges(4,2) * t244;
t204 = t260 * t297 + t240;
t239 = pkin(1) * t263 + pkin(7) * t260;
t280 = -V_base(4) * pkin(6) + t246 * t239 + V_base(2);
t279 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t278 = t241 * t177 + t279;
t277 = pkin(9) * t244 + t245 * t321;
t178 = qJ(3) * t260 + t263 * t322;
t273 = -qJD(3) * t263 + t246 * t178 + t280;
t201 = t289 * t260;
t215 = t244 * pkin(3) - t245 * pkin(8);
t272 = t240 * t215 + (-t201 + t300) * t246 + t290;
t202 = t289 * t263;
t271 = t241 * t201 + (-t178 - t202) * t240 + t278;
t145 = -pkin(4) * t306 + t260 * t277;
t156 = -pkin(9) * t245 + t244 * t321;
t219 = -qJD(4) * t245 + t246;
t270 = -t145 * t219 + t204 * t156 + t272;
t146 = pkin(4) * t307 + t263 * t277;
t269 = t205 * t145 - t146 * t204 + t271;
t268 = t246 * t202 + (-t215 - t324) * t241 + t273;
t267 = t219 * t146 - t156 * t205 + t268;
t251 = Icges(2,4) * t263;
t237 = rSges(2,1) * t263 - rSges(2,2) * t260;
t236 = rSges(2,1) * t260 + rSges(2,2) * t263;
t235 = rSges(3,1) * t259 + rSges(3,2) * t262;
t232 = Icges(2,1) * t263 - t318;
t231 = Icges(2,1) * t260 + t251;
t229 = -Icges(2,2) * t260 + t251;
t228 = Icges(2,2) * t263 + t318;
t222 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t221 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t220 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = rSges(4,1) * t244 + rSges(4,2) * t245;
t209 = t245 * t304 + t307;
t208 = -t245 * t306 + t305;
t207 = t245 * t305 - t306;
t206 = -t245 * t307 - t304;
t203 = (-qJD(4) - qJD(5)) * t245 + t246;
t200 = rSges(3,3) * t260 + t263 * t288;
t199 = -rSges(3,3) * t263 + t260 * t288;
t187 = rSges(4,3) * t260 + t263 * t287;
t186 = -rSges(4,3) * t263 + t260 * t287;
t176 = -rSges(5,3) * t245 + (rSges(5,1) * t261 - rSges(5,2) * t258) * t244;
t175 = V_base(5) * rSges(2,3) - t236 * t246 + t294;
t174 = t237 * t246 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = -Icges(5,5) * t245 + (Icges(5,1) * t261 - Icges(5,4) * t258) * t244;
t172 = -Icges(5,6) * t245 + (Icges(5,4) * t261 - Icges(5,2) * t258) * t244;
t171 = -Icges(5,3) * t245 + (Icges(5,5) * t261 - Icges(5,6) * t258) * t244;
t170 = t263 * t296 + t205;
t169 = t260 * t296 + t204;
t167 = t236 * V_base(4) - t237 * V_base(5) + V_base(3);
t165 = -rSges(6,3) * t245 + (rSges(6,1) * t250 - rSges(6,2) * t249) * t244;
t154 = rSges(5,1) * t209 + rSges(5,2) * t208 + rSges(5,3) * t312;
t153 = rSges(5,1) * t207 + rSges(5,2) * t206 + rSges(5,3) * t313;
t152 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t312;
t151 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t313;
t150 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t312;
t149 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t313;
t148 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t312;
t147 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t313;
t143 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t312;
t141 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t313;
t126 = t235 * t240 + (-t199 - t238) * t246 + t294;
t125 = t200 * t246 - t235 * t241 + t280;
t121 = t199 * t241 - t200 * t240 + t279;
t120 = t214 * t240 + (-t186 + t300) * t246 + t290;
t119 = t187 * t246 + (-t214 - t324) * t241 + t273;
t118 = t186 * t241 + (-t178 - t187) * t240 + t278;
t117 = -t153 * t219 + t176 * t204 + t272;
t116 = t154 * t219 - t176 * t205 + t268;
t115 = t153 * t205 - t154 * t204 + t271;
t114 = -t141 * t203 + t165 * t169 + t270;
t113 = t143 * t203 - t165 * t170 + t267;
t112 = t141 * t170 - t143 * t169 + t269;
t111 = t169 * t301 - t203 * t303 + t263 * t295 + t270;
t110 = -t170 * t301 + t203 * t302 + t260 * t295 + t267;
t109 = -qJD(6) * t245 - t169 * t302 + t170 * t303 + t269;
t1 = m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + t204 * ((t148 * t313 + t150 * t206 + t152 * t207) * t205 + (t147 * t313 + t149 * t206 + t151 * t207) * t204 + (t171 * t313 + t172 * t206 + t173 * t207) * t219) / 0.2e1 + t205 * ((t148 * t312 + t150 * t208 + t152 * t209) * t205 + (t147 * t312 + t149 * t208 + t151 * t209) * t204 + (t171 * t312 + t172 * t208 + t173 * t209) * t219) / 0.2e1 + m(5) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(7) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t219 * ((-t147 * t204 - t148 * t205 - t171 * t219) * t245 + ((-t150 * t258 + t152 * t261) * t205 + (-t149 * t258 + t151 * t261) * t204 + (-t172 * t258 + t173 * t261) * t219) * t244) / 0.2e1 + m(1) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + ((t189 * t334 + t190 * t333 + t313 * t335) * t203 + (t189 * t338 + t190 * t336 + t313 * t340) * t170 + (t339 * t189 + t337 * t190 + t341 * t313) * t169) * t169 / 0.2e1 + ((t191 * t334 + t192 * t333 + t312 * t335) * t203 + (t338 * t191 + t336 * t192 + t340 * t312) * t170 + (t339 * t191 + t337 * t192 + t312 * t341) * t169) * t170 / 0.2e1 + ((-t169 * t341 - t340 * t170 - t335 * t203) * t245 + ((-t249 * t334 + t250 * t333) * t203 + (-t249 * t338 + t250 * t336) * t170 + (-t249 * t339 + t250 * t337) * t169) * t244) * t203 / 0.2e1 + (t332 * t260 - t331 * t263) * t240 / 0.2e1 + (t331 * t260 + t332 * t263) * t241 / 0.2e1 + ((-t228 * t260 + t231 * t263 + Icges(1,4)) * V_base(5) + (-t229 * t260 + t232 * t263 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t228 * t263 + t231 * t260 + Icges(1,2)) * V_base(5) + (t229 * t263 + t232 * t260 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t183 * t245 + t185 * t244 + t196 * t262 + t198 * t259) * t241 + (t182 * t245 + t184 * t244 + t195 * t262 + t197 * t259) * t240 + (t212 * t245 + t213 * t244 + t227 * t262 + t230 * t259 + Icges(2,3)) * t246) * t246 / 0.2e1 + t246 * V_base(4) * (Icges(2,5) * t263 - Icges(2,6) * t260) + t246 * V_base(5) * (Icges(2,5) * t260 + Icges(2,6) * t263) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
