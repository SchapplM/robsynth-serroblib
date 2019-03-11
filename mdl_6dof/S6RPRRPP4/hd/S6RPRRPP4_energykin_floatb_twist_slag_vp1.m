% Calculate kinetic energy for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:36
% EndTime: 2019-03-09 04:38:39
% DurationCPUTime: 3.31s
% Computational Cost: add. (2156->329), mult. (2234->453), div. (0->0), fcn. (2138->10), ass. (0->165)
t340 = Icges(6,1) + Icges(7,1);
t339 = -Icges(6,4) + Icges(7,5);
t338 = Icges(7,4) + Icges(6,5);
t337 = Icges(6,2) + Icges(7,3);
t336 = -Icges(7,6) + Icges(6,6);
t335 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t334 = rSges(7,1) + pkin(5);
t333 = rSges(7,3) + qJ(6);
t255 = pkin(9) + qJ(3);
t248 = cos(t255);
t256 = qJ(4) + pkin(10);
t247 = sin(t256);
t262 = sin(qJ(1));
t305 = t262 * t247;
t249 = cos(t256);
t264 = cos(qJ(1));
t307 = t249 * t264;
t193 = t248 * t305 + t307;
t304 = t262 * t249;
t308 = t247 * t264;
t194 = t248 * t304 - t308;
t246 = sin(t255);
t310 = t246 * t262;
t332 = t337 * t193 + t339 * t194 - t336 * t310;
t195 = t248 * t308 - t304;
t196 = t248 * t307 + t305;
t309 = t246 * t264;
t331 = t337 * t195 + t339 * t196 - t336 * t309;
t330 = t339 * t193 + t340 * t194 + t338 * t310;
t329 = t339 * t195 + t340 * t196 + t338 * t309;
t328 = t336 * t248 + (t337 * t247 + t339 * t249) * t246;
t327 = -t338 * t248 + (t339 * t247 + t340 * t249) * t246;
t263 = cos(qJ(4));
t301 = t263 * t264;
t261 = sin(qJ(4));
t303 = t262 * t261;
t209 = -t248 * t303 - t301;
t302 = t262 * t263;
t306 = t261 * t264;
t210 = t248 * t302 - t306;
t326 = Icges(5,5) * t210 + Icges(5,6) * t209 - t336 * t193 + t338 * t194 - t335 * t310;
t211 = -t248 * t306 + t302;
t212 = t248 * t301 + t303;
t325 = Icges(5,5) * t212 + Icges(5,6) * t211 - t336 * t195 + t338 * t196 - t335 * t309;
t324 = t335 * t248 + (Icges(5,5) * t263 - Icges(5,6) * t261 - t336 * t247 + t338 * t249) * t246;
t257 = sin(pkin(9));
t319 = pkin(2) * t257;
t258 = cos(pkin(9));
t318 = pkin(2) * t258;
t317 = pkin(4) * t263;
t315 = Icges(2,4) * t262;
t314 = Icges(3,4) * t257;
t313 = Icges(3,4) * t258;
t312 = Icges(4,4) * t246;
t311 = Icges(4,4) * t248;
t299 = rSges(7,2) * t310 + t333 * t193 + t194 * t334;
t298 = rSges(7,2) * t309 + t333 * t195 + t196 * t334;
t297 = -rSges(7,2) * t248 + (t333 * t247 + t249 * t334) * t246;
t180 = -pkin(7) * t264 + t262 * t318;
t236 = t262 * pkin(1) - qJ(2) * t264;
t296 = -t180 - t236;
t295 = qJD(4) * t246;
t294 = qJD(5) * t246;
t293 = V_base(4) * t236 + V_base(3);
t292 = V_base(5) * pkin(6) + V_base(1);
t242 = qJD(3) * t262 + V_base(4);
t250 = V_base(6) + qJD(1);
t289 = qJD(2) * t262 + t292;
t288 = V_base(5) * t319 + t289;
t287 = pkin(3) * t248 + pkin(8) * t246;
t241 = -qJD(3) * t264 + V_base(5);
t286 = rSges(3,1) * t258 - rSges(3,2) * t257;
t285 = rSges(4,1) * t248 - rSges(4,2) * t246;
t284 = Icges(3,1) * t258 - t314;
t283 = Icges(4,1) * t248 - t312;
t282 = -Icges(3,2) * t257 + t313;
t281 = -Icges(4,2) * t246 + t311;
t280 = Icges(3,5) * t258 - Icges(3,6) * t257;
t279 = Icges(4,5) * t248 - Icges(4,6) * t246;
t238 = pkin(1) * t264 + t262 * qJ(2);
t278 = -qJD(2) * t264 + t250 * t238 + V_base(2);
t277 = qJ(5) * t246 + t248 * t317;
t276 = (-Icges(4,3) * t264 + t262 * t279) * t241 + (Icges(4,3) * t262 + t264 * t279) * t242 + (Icges(4,5) * t246 + Icges(4,6) * t248) * t250;
t181 = pkin(7) * t262 + t264 * t318;
t275 = V_base(4) * t180 + (-t181 - t238) * V_base(5) + t293;
t274 = (-Icges(3,3) * t264 + t262 * t280) * V_base(5) + (Icges(3,3) * t262 + t264 * t280) * V_base(4) + (Icges(3,5) * t257 + Icges(3,6) * t258) * t250;
t205 = t287 * t262;
t218 = pkin(3) * t246 - pkin(8) * t248;
t273 = t241 * t218 + (-t205 + t296) * t250 + t288;
t206 = t287 * t264;
t272 = t242 * t205 - t206 * t241 + t275;
t271 = t250 * t181 + (-pkin(6) - t319) * V_base(4) + t278;
t161 = -qJ(5) * t248 + t246 * t317;
t207 = t262 * t295 + t241;
t270 = t207 * t161 + t264 * t294 + t273;
t269 = t250 * t206 - t242 * t218 + t271;
t149 = -pkin(4) * t306 + t262 * t277;
t208 = t264 * t295 + t242;
t268 = -qJD(5) * t248 + t208 * t149 + t272;
t150 = pkin(4) * t303 + t264 * t277;
t220 = -qJD(4) * t248 + t250;
t267 = t220 * t150 + t262 * t294 + t269;
t185 = -Icges(4,6) * t264 + t262 * t281;
t186 = Icges(4,6) * t262 + t264 * t281;
t187 = -Icges(4,5) * t264 + t262 * t283;
t188 = Icges(4,5) * t262 + t264 * t283;
t215 = Icges(4,2) * t248 + t312;
t216 = Icges(4,1) * t246 + t311;
t266 = (-t186 * t246 + t188 * t248) * t242 + (-t185 * t246 + t187 * t248) * t241 + (-t215 * t246 + t216 * t248) * t250;
t199 = -Icges(3,6) * t264 + t262 * t282;
t200 = Icges(3,6) * t262 + t264 * t282;
t201 = -Icges(3,5) * t264 + t262 * t284;
t202 = Icges(3,5) * t262 + t264 * t284;
t225 = Icges(3,2) * t258 + t314;
t226 = Icges(3,1) * t257 + t313;
t265 = (-t200 * t257 + t202 * t258) * V_base(4) + (-t199 * t257 + t201 * t258) * V_base(5) + (-t225 * t257 + t226 * t258) * t250;
t253 = Icges(2,4) * t264;
t239 = rSges(2,1) * t264 - t262 * rSges(2,2);
t237 = t262 * rSges(2,1) + rSges(2,2) * t264;
t233 = Icges(2,1) * t264 - t315;
t232 = Icges(2,1) * t262 + t253;
t231 = -Icges(2,2) * t262 + t253;
t230 = Icges(2,2) * t264 + t315;
t229 = Icges(2,5) * t264 - Icges(2,6) * t262;
t228 = Icges(2,5) * t262 + Icges(2,6) * t264;
t227 = rSges(3,1) * t257 + rSges(3,2) * t258;
t223 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t222 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t221 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t217 = rSges(4,1) * t246 + rSges(4,2) * t248;
t204 = t262 * rSges(3,3) + t264 * t286;
t203 = -rSges(3,3) * t264 + t262 * t286;
t191 = t262 * rSges(4,3) + t264 * t285;
t190 = -rSges(4,3) * t264 + t262 * t285;
t179 = -rSges(5,3) * t248 + (rSges(5,1) * t263 - rSges(5,2) * t261) * t246;
t178 = V_base(5) * rSges(2,3) - t237 * t250 + t292;
t177 = t239 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t176 = -Icges(5,5) * t248 + (Icges(5,1) * t263 - Icges(5,4) * t261) * t246;
t175 = -Icges(5,6) * t248 + (Icges(5,4) * t263 - Icges(5,2) * t261) * t246;
t172 = t237 * V_base(4) - t239 * V_base(5) + V_base(3);
t169 = -rSges(6,3) * t248 + (rSges(6,1) * t249 - rSges(6,2) * t247) * t246;
t158 = t212 * rSges(5,1) + t211 * rSges(5,2) + rSges(5,3) * t309;
t157 = rSges(5,1) * t210 + rSges(5,2) * t209 + rSges(5,3) * t310;
t156 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t309;
t155 = Icges(5,1) * t210 + Icges(5,4) * t209 + Icges(5,5) * t310;
t154 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t309;
t153 = Icges(5,4) * t210 + Icges(5,2) * t209 + Icges(5,6) * t310;
t147 = t196 * rSges(6,1) - t195 * rSges(6,2) + rSges(6,3) * t309;
t145 = rSges(6,1) * t194 - rSges(6,2) * t193 + rSges(6,3) * t310;
t130 = t227 * V_base(5) + (-t203 - t236) * t250 + t289;
t129 = t250 * t204 + (-pkin(6) - t227) * V_base(4) + t278;
t127 = t203 * V_base(4) + (-t204 - t238) * V_base(5) + t293;
t126 = t217 * t241 + (-t190 + t296) * t250 + t288;
t125 = t250 * t191 - t242 * t217 + t271;
t124 = t190 * t242 - t191 * t241 + t275;
t123 = -t157 * t220 + t179 * t207 + t273;
t122 = t220 * t158 - t208 * t179 + t269;
t121 = t157 * t208 - t158 * t207 + t272;
t120 = t169 * t207 + (-t145 - t149) * t220 + t270;
t119 = t220 * t147 + (-t161 - t169) * t208 + t267;
t118 = t145 * t208 + (-t147 - t150) * t207 + t268;
t117 = qJD(6) * t195 + t297 * t207 + (-t149 - t299) * t220 + t270;
t116 = qJD(6) * t193 + t298 * t220 + (-t161 - t297) * t208 + t267;
t115 = qJD(6) * t246 * t247 + t299 * t208 + (-t150 - t298) * t207 + t268;
t1 = t242 * (t276 * t262 + t266 * t264) / 0.2e1 + t241 * (t266 * t262 - t276 * t264) / 0.2e1 + m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(2) * (t172 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(4) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(6) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(7) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + ((t175 * t209 + t176 * t210 + t193 * t328 + t194 * t327 + t310 * t324) * t220 + (t154 * t209 + t156 * t210 + t193 * t331 + t194 * t329 + t310 * t325) * t208 + (t209 * t153 + t210 * t155 + t332 * t193 + t330 * t194 + t326 * t310) * t207) * t207 / 0.2e1 + ((t211 * t175 + t212 * t176 + t195 * t328 + t196 * t327 + t309 * t324) * t220 + (t211 * t154 + t212 * t156 + t331 * t195 + t329 * t196 + t325 * t309) * t208 + (t211 * t153 + t212 * t155 + t195 * t332 + t330 * t196 + t326 * t309) * t207) * t208 / 0.2e1 + ((-t207 * t326 - t208 * t325 - t220 * t324) * t248 + ((-t175 * t261 + t176 * t263 + t247 * t328 + t249 * t327) * t220 + (-t154 * t261 + t156 * t263 + t247 * t331 + t249 * t329) * t208 + (-t153 * t261 + t155 * t263 + t247 * t332 + t330 * t249) * t207) * t246) * t220 / 0.2e1 + ((t186 * t248 + t188 * t246) * t242 + (t185 * t248 + t187 * t246) * t241 + (t199 * t258 + t201 * t257 + t228) * V_base(5) + (t200 * t258 + t202 * t257 + t229) * V_base(4) + (t248 * t215 + t246 * t216 + t258 * t225 + t257 * t226 + Icges(2,3)) * t250) * t250 / 0.2e1 + (t229 * t250 + t262 * t274 + t264 * t265 + (-t262 * t230 + t232 * t264 + Icges(1,4)) * V_base(5) + (-t262 * t231 + t264 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t228 * t250 + t262 * t265 - t264 * t274 + (t264 * t230 + t262 * t232 + Icges(1,2)) * V_base(5) + (t231 * t264 + t262 * t233 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
