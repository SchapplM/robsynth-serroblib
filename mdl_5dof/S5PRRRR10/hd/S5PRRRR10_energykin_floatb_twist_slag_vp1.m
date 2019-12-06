% Calculate kinetic energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:10
% EndTime: 2019-12-05 17:23:13
% DurationCPUTime: 3.40s
% Computational Cost: add. (3562->365), mult. (9334->553), div. (0->0), fcn. (11822->14), ass. (0->160)
t346 = cos(qJ(3));
t345 = cos(qJ(4));
t309 = cos(pkin(5));
t344 = pkin(7) * t309;
t343 = cos(pkin(6));
t342 = sin(pkin(6));
t306 = sin(pkin(11));
t341 = Icges(2,4) * t306;
t307 = sin(pkin(5));
t340 = t306 * t307;
t308 = cos(pkin(11));
t339 = t307 * t308;
t313 = sin(qJ(2));
t338 = t307 * t313;
t337 = t309 * t313;
t315 = cos(qJ(2));
t336 = t309 * t315;
t335 = qJD(2) * t307;
t334 = V_base(5) * qJ(1) + V_base(1);
t330 = qJD(1) + V_base(3);
t290 = t306 * t335 + V_base(4);
t300 = qJD(2) * t309 + V_base(6);
t329 = t307 * t343;
t328 = t307 * t342;
t278 = -t306 * t336 - t308 * t313;
t262 = -t278 * t342 + t306 * t329;
t251 = qJD(3) * t262 + t290;
t275 = t309 * t343 - t315 * t328;
t263 = qJD(3) * t275 + t300;
t327 = t343 * t346;
t326 = t346 * t342;
t279 = -t306 * t337 + t308 * t315;
t312 = sin(qJ(3));
t325 = t307 * t326;
t236 = -t278 * t327 + t279 * t312 - t306 * t325;
t217 = qJD(4) * t236 + t251;
t259 = -t307 * t315 * t327 - t309 * t326 + t312 * t338;
t230 = qJD(4) * t259 + t263;
t289 = -t308 * t335 + V_base(5);
t276 = -t306 * t313 + t308 * t336;
t261 = -t276 * t342 - t308 * t329;
t250 = qJD(3) * t261 + t289;
t282 = pkin(1) * t306 - pkin(7) * t339;
t324 = -t282 * V_base(6) + t344 * V_base(5) + t334;
t283 = pkin(1) * t308 + pkin(7) * t340;
t323 = t282 * V_base(4) - t283 * V_base(5) + t330;
t277 = t306 * t315 + t308 * t337;
t234 = -t276 * t327 + t277 * t312 + t308 * t325;
t216 = qJD(4) * t234 + t250;
t322 = V_base(6) * t283 + V_base(2) + (-qJ(1) - t344) * V_base(4);
t240 = pkin(2) * t277 + pkin(8) * t261;
t264 = pkin(2) * t338 + pkin(8) * t275;
t321 = -t240 * t300 + t264 * t289 + t324;
t241 = pkin(2) * t279 + pkin(8) * t262;
t320 = t240 * t290 - t241 * t289 + t323;
t319 = t241 * t300 - t264 * t290 + t322;
t235 = t277 * t346 + (t276 * t343 - t308 * t328) * t312;
t212 = pkin(3) * t235 + pkin(9) * t234;
t260 = t309 * t342 * t312 + (t312 * t315 * t343 + t313 * t346) * t307;
t228 = pkin(3) * t260 + pkin(9) * t259;
t318 = -t212 * t263 + t228 * t250 + t321;
t237 = t279 * t346 + (t278 * t343 + t306 * t328) * t312;
t213 = pkin(3) * t237 + pkin(9) * t236;
t317 = t212 * t251 - t213 * t250 + t320;
t316 = t213 * t263 - t228 * t251 + t319;
t314 = cos(qJ(5));
t311 = sin(qJ(4));
t310 = sin(qJ(5));
t304 = Icges(2,4) * t308;
t298 = rSges(2,1) * t308 - rSges(2,2) * t306;
t297 = rSges(2,1) * t306 + rSges(2,2) * t308;
t296 = Icges(2,1) * t308 - t341;
t295 = Icges(2,1) * t306 + t304;
t294 = -Icges(2,2) * t306 + t304;
t293 = Icges(2,2) * t308 + t341;
t288 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t287 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t286 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t272 = t309 * rSges(3,3) + (rSges(3,1) * t313 + rSges(3,2) * t315) * t307;
t271 = Icges(3,5) * t309 + (Icges(3,1) * t313 + Icges(3,4) * t315) * t307;
t270 = Icges(3,6) * t309 + (Icges(3,4) * t313 + Icges(3,2) * t315) * t307;
t269 = Icges(3,3) * t309 + (Icges(3,5) * t313 + Icges(3,6) * t315) * t307;
t266 = V_base(5) * rSges(2,3) - t297 * V_base(6) + t334;
t265 = t298 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t258 = t297 * V_base(4) - t298 * V_base(5) + t330;
t249 = rSges(3,1) * t279 + rSges(3,2) * t278 + rSges(3,3) * t340;
t248 = rSges(3,1) * t277 + rSges(3,2) * t276 - rSges(3,3) * t339;
t247 = Icges(3,1) * t279 + Icges(3,4) * t278 + Icges(3,5) * t340;
t246 = Icges(3,1) * t277 + Icges(3,4) * t276 - Icges(3,5) * t339;
t245 = Icges(3,4) * t279 + Icges(3,2) * t278 + Icges(3,6) * t340;
t244 = Icges(3,4) * t277 + Icges(3,2) * t276 - Icges(3,6) * t339;
t243 = Icges(3,5) * t279 + Icges(3,6) * t278 + Icges(3,3) * t340;
t242 = Icges(3,5) * t277 + Icges(3,6) * t276 - Icges(3,3) * t339;
t239 = t260 * t345 + t275 * t311;
t238 = t260 * t311 - t275 * t345;
t227 = rSges(4,1) * t260 - rSges(4,2) * t259 + rSges(4,3) * t275;
t226 = Icges(4,1) * t260 - Icges(4,4) * t259 + Icges(4,5) * t275;
t225 = Icges(4,4) * t260 - Icges(4,2) * t259 + Icges(4,6) * t275;
t224 = Icges(4,5) * t260 - Icges(4,6) * t259 + Icges(4,3) * t275;
t223 = t237 * t345 + t262 * t311;
t222 = t237 * t311 - t262 * t345;
t221 = t235 * t345 + t261 * t311;
t220 = t235 * t311 - t261 * t345;
t219 = t239 * t314 + t259 * t310;
t218 = -t239 * t310 + t259 * t314;
t214 = pkin(4) * t239 + pkin(10) * t238;
t211 = -t248 * t300 + t272 * t289 + t324;
t210 = t249 * t300 - t272 * t290 + t322;
t209 = qJD(5) * t238 + t230;
t207 = rSges(5,1) * t239 - rSges(5,2) * t238 + rSges(5,3) * t259;
t206 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t262;
t205 = rSges(4,1) * t235 - rSges(4,2) * t234 + rSges(4,3) * t261;
t204 = Icges(5,1) * t239 - Icges(5,4) * t238 + Icges(5,5) * t259;
t203 = Icges(5,4) * t239 - Icges(5,2) * t238 + Icges(5,6) * t259;
t202 = Icges(5,5) * t239 - Icges(5,6) * t238 + Icges(5,3) * t259;
t201 = Icges(4,1) * t237 - Icges(4,4) * t236 + Icges(4,5) * t262;
t200 = Icges(4,1) * t235 - Icges(4,4) * t234 + Icges(4,5) * t261;
t199 = Icges(4,4) * t237 - Icges(4,2) * t236 + Icges(4,6) * t262;
t198 = Icges(4,4) * t235 - Icges(4,2) * t234 + Icges(4,6) * t261;
t197 = Icges(4,5) * t237 - Icges(4,6) * t236 + Icges(4,3) * t262;
t196 = Icges(4,5) * t235 - Icges(4,6) * t234 + Icges(4,3) * t261;
t195 = t248 * t290 - t249 * t289 + t323;
t194 = t223 * t314 + t236 * t310;
t193 = -t223 * t310 + t236 * t314;
t192 = t221 * t314 + t234 * t310;
t191 = -t221 * t310 + t234 * t314;
t189 = pkin(4) * t223 + pkin(10) * t222;
t188 = pkin(4) * t221 + pkin(10) * t220;
t187 = qJD(5) * t222 + t217;
t186 = qJD(5) * t220 + t216;
t185 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t238;
t184 = rSges(5,1) * t223 - rSges(5,2) * t222 + rSges(5,3) * t236;
t183 = rSges(5,1) * t221 - rSges(5,2) * t220 + rSges(5,3) * t234;
t182 = Icges(5,1) * t223 - Icges(5,4) * t222 + Icges(5,5) * t236;
t181 = Icges(5,1) * t221 - Icges(5,4) * t220 + Icges(5,5) * t234;
t180 = Icges(6,1) * t219 + Icges(6,4) * t218 + Icges(6,5) * t238;
t179 = Icges(5,4) * t223 - Icges(5,2) * t222 + Icges(5,6) * t236;
t178 = Icges(5,4) * t221 - Icges(5,2) * t220 + Icges(5,6) * t234;
t177 = Icges(6,4) * t219 + Icges(6,2) * t218 + Icges(6,6) * t238;
t176 = Icges(5,5) * t223 - Icges(5,6) * t222 + Icges(5,3) * t236;
t175 = Icges(5,5) * t221 - Icges(5,6) * t220 + Icges(5,3) * t234;
t174 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t238;
t173 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t222;
t172 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t220;
t171 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t222;
t170 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t220;
t169 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t222;
t168 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t220;
t167 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t222;
t166 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t220;
t165 = -t205 * t263 + t227 * t250 + t321;
t164 = t206 * t263 - t227 * t251 + t319;
t163 = t205 * t251 - t206 * t250 + t320;
t162 = -t183 * t230 + t207 * t216 + t318;
t161 = t184 * t230 - t207 * t217 + t316;
t160 = t183 * t217 - t184 * t216 + t317;
t159 = -t172 * t209 + t185 * t186 - t188 * t230 + t214 * t216 + t318;
t158 = t173 * t209 - t185 * t187 + t189 * t230 - t214 * t217 + t316;
t157 = t172 * t187 - t173 * t186 + t188 * t217 - t189 * t216 + t317;
t1 = m(1) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(2) * (t258 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(3) * (t195 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + t290 * ((t243 * t340 + t245 * t278 + t247 * t279) * t290 + (t242 * t340 + t244 * t278 + t246 * t279) * t289 + (t269 * t340 + t270 * t278 + t271 * t279) * t300) / 0.2e1 + t289 * ((-t243 * t339 + t245 * t276 + t247 * t277) * t290 + (-t242 * t339 + t244 * t276 + t246 * t277) * t289 + (-t269 * t339 + t270 * t276 + t271 * t277) * t300) / 0.2e1 + t300 * ((t242 * t289 + t243 * t290 + t269 * t300) * t309 + ((t245 * t315 + t247 * t313) * t290 + (t244 * t315 + t246 * t313) * t289 + (t270 * t315 + t271 * t313) * t300) * t307) / 0.2e1 + m(4) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + t251 * ((t262 * t197 - t236 * t199 + t237 * t201) * t251 + (t196 * t262 - t198 * t236 + t200 * t237) * t250 + (t224 * t262 - t225 * t236 + t226 * t237) * t263) / 0.2e1 + t250 * ((t197 * t261 - t199 * t234 + t201 * t235) * t251 + (t261 * t196 - t234 * t198 + t235 * t200) * t250 + (t224 * t261 - t225 * t234 + t226 * t235) * t263) / 0.2e1 + t263 * ((t197 * t275 - t199 * t259 + t201 * t260) * t251 + (t196 * t275 - t198 * t259 + t200 * t260) * t250 + (t224 * t275 - t225 * t259 + t226 * t260) * t263) / 0.2e1 + m(5) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + t217 * ((t236 * t176 - t222 * t179 + t223 * t182) * t217 + (t175 * t236 - t178 * t222 + t181 * t223) * t216 + (t202 * t236 - t203 * t222 + t204 * t223) * t230) / 0.2e1 + t216 * ((t176 * t234 - t179 * t220 + t182 * t221) * t217 + (t234 * t175 - t220 * t178 + t221 * t181) * t216 + (t202 * t234 - t203 * t220 + t204 * t221) * t230) / 0.2e1 + t230 * ((t176 * t259 - t179 * t238 + t182 * t239) * t217 + (t175 * t259 - t178 * t238 + t181 * t239) * t216 + (t259 * t202 - t238 * t203 + t239 * t204) * t230) / 0.2e1 + m(6) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + t187 * ((t222 * t167 + t193 * t169 + t194 * t171) * t187 + (t166 * t222 + t168 * t193 + t170 * t194) * t186 + (t174 * t222 + t177 * t193 + t180 * t194) * t209) / 0.2e1 + t186 * ((t167 * t220 + t169 * t191 + t171 * t192) * t187 + (t220 * t166 + t191 * t168 + t192 * t170) * t186 + (t174 * t220 + t177 * t191 + t180 * t192) * t209) / 0.2e1 + t209 * ((t167 * t238 + t169 * t218 + t171 * t219) * t187 + (t166 * t238 + t168 * t218 + t170 * t219) * t186 + (t238 * t174 + t218 * t177 + t219 * t180) * t209) / 0.2e1 + ((-t293 * t306 + t295 * t308 + Icges(1,4)) * V_base(5) + (-t294 * t306 + t296 * t308 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t293 * t308 + t295 * t306 + Icges(1,2)) * V_base(5) + (t294 * t308 + t296 * t306 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t306 + Icges(2,6) * t308 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t308 - Icges(2,6) * t306 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
