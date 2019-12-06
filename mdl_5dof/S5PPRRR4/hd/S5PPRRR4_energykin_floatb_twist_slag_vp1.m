% Calculate kinetic energy for
% S5PPRRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:32
% EndTime: 2019-12-05 15:18:35
% DurationCPUTime: 2.94s
% Computational Cost: add. (3466->370), mult. (9174->537), div. (0->0), fcn. (11662->14), ass. (0->160)
t304 = sin(pkin(10));
t307 = cos(pkin(10));
t346 = Icges(2,5) * t307 - Icges(2,6) * t304 + Icges(1,5);
t345 = Icges(2,5) * t304 + Icges(2,6) * t307 + Icges(1,6);
t344 = cos(qJ(3));
t343 = cos(qJ(4));
t342 = cos(pkin(6));
t341 = sin(pkin(6));
t340 = Icges(2,4) * t304;
t308 = cos(pkin(5));
t339 = qJ(2) * t308;
t303 = sin(pkin(11));
t305 = sin(pkin(5));
t338 = t303 * t305;
t337 = t304 * t305;
t336 = t304 * t308;
t335 = t305 * t307;
t334 = t307 * t308;
t333 = qJD(2) * t305;
t332 = V_base(5) * qJ(1) + V_base(1);
t328 = qJD(1) + V_base(3);
t306 = cos(pkin(11));
t278 = -t303 * t307 - t306 * t336;
t327 = t305 * t342;
t261 = -t278 * t341 + t304 * t327;
t251 = qJD(3) * t261 + V_base(4);
t276 = -t303 * t304 + t306 * t334;
t260 = -t276 * t341 - t307 * t327;
t250 = qJD(3) * t260 + V_base(5);
t326 = t305 * t341;
t275 = -t306 * t326 + t308 * t342;
t268 = qJD(3) * t275 + V_base(6);
t325 = -qJ(1) - t339;
t279 = -t303 * t336 + t306 * t307;
t311 = sin(qJ(3));
t323 = t344 * t341;
t320 = t305 * t323;
t324 = t342 * t344;
t233 = -t278 * t324 + t279 * t311 - t304 * t320;
t217 = qJD(4) * t233 + t251;
t277 = t303 * t334 + t304 * t306;
t231 = -t276 * t324 + t277 * t311 + t307 * t320;
t216 = qJD(4) * t231 + t250;
t257 = -t305 * t306 * t324 - t308 * t323 + t311 * t338;
t235 = qJD(4) * t257 + t268;
t322 = t304 * t333 + V_base(5) * t339 + t332;
t282 = pkin(1) * t304 - qJ(2) * t335;
t321 = qJD(2) * t308 + V_base(4) * t282 + t328;
t283 = pkin(1) * t307 + qJ(2) * t337;
t319 = V_base(6) * t283 - t307 * t333 + V_base(2);
t240 = t277 * pkin(2) + pkin(7) * t260;
t263 = pkin(2) * t338 + pkin(7) * t275;
t318 = V_base(5) * t263 + (-t240 - t282) * V_base(6) + t322;
t241 = t279 * pkin(2) + pkin(7) * t261;
t317 = V_base(4) * t240 + (-t241 - t283) * V_base(5) + t321;
t232 = t277 * t344 + (t276 * t342 - t307 * t326) * t311;
t209 = pkin(3) * t232 + pkin(8) * t231;
t258 = t308 * t341 * t311 + (t342 * t306 * t311 + t303 * t344) * t305;
t228 = pkin(3) * t258 + pkin(8) * t257;
t316 = -t209 * t268 + t250 * t228 + t318;
t234 = t279 * t344 + (t278 * t342 + t304 * t326) * t311;
t210 = pkin(3) * t234 + pkin(8) * t233;
t315 = t251 * t209 - t210 * t250 + t317;
t314 = V_base(6) * t241 + (-t263 + t325) * V_base(4) + t319;
t313 = t268 * t210 - t228 * t251 + t314;
t312 = cos(qJ(5));
t310 = sin(qJ(4));
t309 = sin(qJ(5));
t301 = Icges(2,4) * t307;
t296 = rSges(2,1) * t307 - rSges(2,2) * t304;
t295 = rSges(2,1) * t304 + rSges(2,2) * t307;
t294 = Icges(2,1) * t307 - t340;
t293 = Icges(2,1) * t304 + t301;
t292 = -Icges(2,2) * t304 + t301;
t291 = Icges(2,2) * t307 + t340;
t288 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t287 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t286 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t272 = rSges(3,3) * t308 + (rSges(3,1) * t303 + rSges(3,2) * t306) * t305;
t271 = Icges(3,5) * t308 + (Icges(3,1) * t303 + Icges(3,4) * t306) * t305;
t270 = Icges(3,6) * t308 + (Icges(3,4) * t303 + Icges(3,2) * t306) * t305;
t269 = Icges(3,3) * t308 + (Icges(3,5) * t303 + Icges(3,6) * t306) * t305;
t265 = V_base(5) * rSges(2,3) - t295 * V_base(6) + t332;
t264 = t296 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t259 = t295 * V_base(4) - t296 * V_base(5) + t328;
t249 = rSges(3,1) * t279 + rSges(3,2) * t278 + rSges(3,3) * t337;
t248 = rSges(3,1) * t277 + rSges(3,2) * t276 - rSges(3,3) * t335;
t247 = Icges(3,1) * t279 + Icges(3,4) * t278 + Icges(3,5) * t337;
t246 = Icges(3,1) * t277 + Icges(3,4) * t276 - Icges(3,5) * t335;
t245 = Icges(3,4) * t279 + Icges(3,2) * t278 + Icges(3,6) * t337;
t244 = Icges(3,4) * t277 + Icges(3,2) * t276 - Icges(3,6) * t335;
t243 = Icges(3,5) * t279 + Icges(3,6) * t278 + Icges(3,3) * t337;
t242 = Icges(3,5) * t277 + Icges(3,6) * t276 - Icges(3,3) * t335;
t237 = t258 * t343 + t275 * t310;
t236 = t258 * t310 - t275 * t343;
t227 = rSges(4,1) * t258 - rSges(4,2) * t257 + rSges(4,3) * t275;
t226 = Icges(4,1) * t258 - Icges(4,4) * t257 + Icges(4,5) * t275;
t225 = Icges(4,4) * t258 - Icges(4,2) * t257 + Icges(4,6) * t275;
t224 = Icges(4,5) * t258 - Icges(4,6) * t257 + Icges(4,3) * t275;
t223 = t237 * t312 + t257 * t309;
t222 = -t237 * t309 + t257 * t312;
t221 = t234 * t343 + t261 * t310;
t220 = t234 * t310 - t261 * t343;
t219 = t232 * t343 + t260 * t310;
t218 = t232 * t310 - t260 * t343;
t214 = t272 * V_base(5) + (-t248 - t282) * V_base(6) + t322;
t213 = t249 * V_base(6) + (-t272 + t325) * V_base(4) + t319;
t212 = qJD(5) * t236 + t235;
t211 = pkin(4) * t237 + pkin(9) * t236;
t207 = t248 * V_base(4) + (-t249 - t283) * V_base(5) + t321;
t206 = rSges(5,1) * t237 - rSges(5,2) * t236 + rSges(5,3) * t257;
t205 = Icges(5,1) * t237 - Icges(5,4) * t236 + Icges(5,5) * t257;
t204 = Icges(5,4) * t237 - Icges(5,2) * t236 + Icges(5,6) * t257;
t203 = Icges(5,5) * t237 - Icges(5,6) * t236 + Icges(5,3) * t257;
t202 = rSges(4,1) * t234 - rSges(4,2) * t233 + rSges(4,3) * t261;
t201 = rSges(4,1) * t232 - rSges(4,2) * t231 + rSges(4,3) * t260;
t200 = Icges(4,1) * t234 - Icges(4,4) * t233 + Icges(4,5) * t261;
t199 = Icges(4,1) * t232 - Icges(4,4) * t231 + Icges(4,5) * t260;
t198 = Icges(4,4) * t234 - Icges(4,2) * t233 + Icges(4,6) * t261;
t197 = Icges(4,4) * t232 - Icges(4,2) * t231 + Icges(4,6) * t260;
t196 = Icges(4,5) * t234 - Icges(4,6) * t233 + Icges(4,3) * t261;
t195 = Icges(4,5) * t232 - Icges(4,6) * t231 + Icges(4,3) * t260;
t193 = t221 * t312 + t233 * t309;
t192 = -t221 * t309 + t233 * t312;
t191 = t219 * t312 + t231 * t309;
t190 = -t219 * t309 + t231 * t312;
t189 = qJD(5) * t220 + t217;
t188 = qJD(5) * t218 + t216;
t187 = pkin(4) * t221 + pkin(9) * t220;
t186 = pkin(4) * t219 + pkin(9) * t218;
t185 = rSges(6,1) * t223 + rSges(6,2) * t222 + rSges(6,3) * t236;
t184 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t236;
t183 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t236;
t182 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t236;
t181 = rSges(5,1) * t221 - rSges(5,2) * t220 + rSges(5,3) * t233;
t180 = rSges(5,1) * t219 - rSges(5,2) * t218 + rSges(5,3) * t231;
t179 = Icges(5,1) * t221 - Icges(5,4) * t220 + Icges(5,5) * t233;
t178 = Icges(5,1) * t219 - Icges(5,4) * t218 + Icges(5,5) * t231;
t177 = Icges(5,4) * t221 - Icges(5,2) * t220 + Icges(5,6) * t233;
t176 = Icges(5,4) * t219 - Icges(5,2) * t218 + Icges(5,6) * t231;
t175 = Icges(5,5) * t221 - Icges(5,6) * t220 + Icges(5,3) * t233;
t174 = Icges(5,5) * t219 - Icges(5,6) * t218 + Icges(5,3) * t231;
t173 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t220;
t172 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t218;
t171 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t220;
t170 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t218;
t169 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t220;
t168 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t218;
t167 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t220;
t166 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t218;
t165 = -t201 * t268 + t227 * t250 + t318;
t164 = t202 * t268 - t227 * t251 + t314;
t163 = t201 * t251 - t202 * t250 + t317;
t162 = -t180 * t235 + t206 * t216 + t316;
t161 = t181 * t235 - t206 * t217 + t313;
t160 = t180 * t217 - t181 * t216 + t315;
t159 = -t172 * t212 + t185 * t188 - t186 * t235 + t211 * t216 + t316;
t158 = t173 * t212 - t185 * t189 + t187 * t235 - t211 * t217 + t313;
t157 = t172 * t189 - t173 * t188 + t186 * t217 - t187 * t216 + t315;
t1 = m(1) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + m(2) * (t259 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(3) * (t207 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + m(4) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + t251 * ((t261 * t196 - t233 * t198 + t234 * t200) * t251 + (t195 * t261 - t197 * t233 + t199 * t234) * t250 + (t224 * t261 - t225 * t233 + t226 * t234) * t268) / 0.2e1 + t250 * ((t196 * t260 - t198 * t231 + t200 * t232) * t251 + (t260 * t195 - t231 * t197 + t232 * t199) * t250 + (t224 * t260 - t225 * t231 + t226 * t232) * t268) / 0.2e1 + t268 * ((t196 * t275 - t198 * t257 + t200 * t258) * t251 + (t195 * t275 - t197 * t257 + t199 * t258) * t250 + (t224 * t275 - t225 * t257 + t226 * t258) * t268) / 0.2e1 + m(5) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + t217 * ((t233 * t175 - t220 * t177 + t221 * t179) * t217 + (t174 * t233 - t176 * t220 + t178 * t221) * t216 + (t203 * t233 - t204 * t220 + t205 * t221) * t235) / 0.2e1 + t216 * ((t175 * t231 - t177 * t218 + t179 * t219) * t217 + (t231 * t174 - t218 * t176 + t219 * t178) * t216 + (t203 * t231 - t204 * t218 + t205 * t219) * t235) / 0.2e1 + t235 * ((t175 * t257 - t177 * t236 + t179 * t237) * t217 + (t174 * t257 - t176 * t236 + t178 * t237) * t216 + (t257 * t203 - t236 * t204 + t237 * t205) * t235) / 0.2e1 + m(6) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + t189 * ((t220 * t167 + t192 * t169 + t193 * t171) * t189 + (t166 * t220 + t168 * t192 + t170 * t193) * t188 + (t182 * t220 + t183 * t192 + t184 * t193) * t212) / 0.2e1 + t188 * ((t167 * t218 + t169 * t190 + t171 * t191) * t189 + (t218 * t166 + t190 * t168 + t191 * t170) * t188 + (t182 * t218 + t183 * t190 + t184 * t191) * t212) / 0.2e1 + t212 * ((t167 * t236 + t169 * t222 + t171 * t223) * t189 + (t166 * t236 + t168 * t222 + t170 * t223) * t188 + (t236 * t182 + t222 * t183 + t223 * t184) * t212) / 0.2e1 + ((t269 * t337 + t270 * t278 + t271 * t279 + t346) * V_base(6) + (t242 * t337 + t244 * t278 + t246 * t279 - t291 * t304 + t293 * t307 + Icges(1,4)) * V_base(5) + (t243 * t337 + t245 * t278 + t247 * t279 - t292 * t304 + t294 * t307 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t269 * t335 + t270 * t276 + t271 * t277 + t345) * V_base(6) + (-t242 * t335 + t244 * t276 + t246 * t277 + t291 * t307 + t293 * t304 + Icges(1,2)) * V_base(5) + (-t243 * t335 + t245 * t276 + t247 * t277 + t292 * t307 + t294 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(1,3) + Icges(2,3) + t269 * t308 + (t270 * t306 + t271 * t303) * t305) * V_base(6) + (t242 * t308 + (t244 * t306 + t246 * t303) * t305 + t345) * V_base(5) + (t243 * t308 + (t245 * t306 + t247 * t303) * t305 + t346) * V_base(4)) * V_base(6) / 0.2e1;
T = t1;
