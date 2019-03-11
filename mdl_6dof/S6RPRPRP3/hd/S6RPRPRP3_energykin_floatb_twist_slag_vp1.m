% Calculate kinetic energy for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:44
% EndTime: 2019-03-09 03:07:46
% DurationCPUTime: 2.95s
% Computational Cost: add. (2124->319), mult. (1989->443), div. (0->0), fcn. (1891->10), ass. (0->154)
t320 = Icges(6,1) + Icges(7,1);
t319 = -Icges(6,4) + Icges(7,5);
t318 = Icges(7,4) + Icges(6,5);
t317 = Icges(6,2) + Icges(7,3);
t316 = -Icges(7,6) + Icges(6,6);
t315 = -Icges(6,3) - Icges(7,2);
t314 = rSges(7,1) + pkin(5);
t313 = rSges(7,3) + qJ(6);
t241 = pkin(10) + qJ(5);
t233 = sin(t241);
t235 = cos(t241);
t242 = qJ(1) + pkin(9);
t236 = cos(t242);
t234 = sin(t242);
t248 = cos(qJ(3));
t285 = t234 * t248;
t168 = t233 * t285 + t235 * t236;
t169 = -t233 * t236 + t235 * t285;
t246 = sin(qJ(3));
t286 = t234 * t246;
t312 = t168 * t317 + t169 * t319 - t286 * t316;
t282 = t236 * t248;
t170 = t233 * t282 - t234 * t235;
t171 = t233 * t234 + t235 * t282;
t283 = t236 * t246;
t311 = t170 * t317 + t171 * t319 - t283 * t316;
t310 = -t168 * t316 + t169 * t318 - t286 * t315;
t309 = -t170 * t316 + t171 * t318 - t283 * t315;
t308 = t319 * t168 + t169 * t320 + t318 * t286;
t307 = t319 * t170 + t171 * t320 + t318 * t283;
t306 = t316 * t248 + (t233 * t317 + t235 * t319) * t246;
t305 = t315 * t248 + (-t233 * t316 + t235 * t318) * t246;
t304 = -t318 * t248 + (t319 * t233 + t235 * t320) * t246;
t247 = sin(qJ(1));
t295 = pkin(1) * t247;
t249 = cos(qJ(1));
t294 = pkin(1) * t249;
t244 = cos(pkin(10));
t293 = pkin(4) * t244;
t292 = -pkin(6) - qJ(2);
t291 = Icges(2,4) * t247;
t290 = Icges(3,4) * t234;
t289 = Icges(4,4) * t246;
t288 = Icges(4,4) * t248;
t243 = sin(pkin(10));
t287 = t234 * t243;
t284 = t236 * t243;
t281 = t243 * t248;
t280 = t244 * t248;
t278 = rSges(7,2) * t286 + t313 * t168 + t314 * t169;
t277 = rSges(7,2) * t283 + t313 * t170 + t314 * t171;
t276 = -rSges(7,2) * t248 + (t313 * t233 + t314 * t235) * t246;
t275 = qJD(4) * t246;
t274 = qJD(5) * t246;
t237 = V_base(6) + qJD(1);
t273 = t237 * t294 + V_base(2);
t272 = V_base(5) * pkin(6) + V_base(1);
t211 = qJD(3) * t234 + V_base(4);
t204 = pkin(2) * t234 - pkin(7) * t236;
t269 = -t204 - t295;
t268 = V_base(5) * qJ(2) + t272;
t267 = V_base(4) * t295 + qJD(2) + V_base(3);
t264 = pkin(3) * t248 + qJ(4) * t246;
t191 = t264 * t234;
t266 = -t191 + t269;
t210 = -qJD(3) * t236 + V_base(5);
t265 = rSges(4,1) * t248 - rSges(4,2) * t246;
t263 = Icges(4,1) * t248 - t289;
t262 = -Icges(4,2) * t246 + t288;
t261 = Icges(4,5) * t248 - Icges(4,6) * t246;
t224 = pkin(3) * t246 - qJ(4) * t248;
t260 = t210 * t224 + t236 * t275 + t268;
t259 = (-Icges(4,3) * t236 + t234 * t261) * t210 + (Icges(4,3) * t234 + t236 * t261) * t211 + (Icges(4,5) * t246 + Icges(4,6) * t248) * t237;
t258 = pkin(8) * t246 + t248 * t293;
t205 = pkin(2) * t236 + pkin(7) * t234;
t257 = t237 * t205 + t292 * V_base(4) + t273;
t193 = t264 * t236;
t256 = t237 * t193 + t234 * t275 + t257;
t255 = V_base(4) * t204 + (-t205 - t294) * V_base(5) + t267;
t146 = -pkin(4) * t284 + t234 * t258;
t166 = -pkin(8) * t248 + t246 * t293;
t254 = t210 * t166 + (-t146 + t266) * t237 + t260;
t253 = -qJD(4) * t248 + t211 * t191 + t255;
t147 = pkin(4) * t287 + t236 * t258;
t252 = t237 * t147 + (-t166 - t224) * t211 + t256;
t251 = t211 * t146 + (-t147 - t193) * t210 + t253;
t158 = -Icges(4,6) * t236 + t234 * t262;
t159 = Icges(4,6) * t234 + t236 * t262;
t160 = -Icges(4,5) * t236 + t234 * t263;
t161 = Icges(4,5) * t234 + t236 * t263;
t215 = Icges(4,2) * t248 + t289;
t218 = Icges(4,1) * t246 + t288;
t250 = (-t159 * t246 + t161 * t248) * t211 + (-t158 * t246 + t160 * t248) * t210 + (-t215 * t246 + t218 * t248) * t237;
t239 = Icges(2,4) * t249;
t231 = Icges(3,4) * t236;
t227 = rSges(2,1) * t249 - t247 * rSges(2,2);
t226 = t247 * rSges(2,1) + rSges(2,2) * t249;
t225 = rSges(4,1) * t246 + rSges(4,2) * t248;
t223 = -qJD(5) * t248 + t237;
t220 = Icges(2,1) * t249 - t291;
t219 = Icges(2,1) * t247 + t239;
t217 = -Icges(2,2) * t247 + t239;
t216 = Icges(2,2) * t249 + t291;
t208 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t207 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t206 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t203 = rSges(3,1) * t236 - rSges(3,2) * t234;
t202 = rSges(3,1) * t234 + rSges(3,2) * t236;
t201 = Icges(3,1) * t236 - t290;
t200 = Icges(3,1) * t234 + t231;
t199 = -Icges(3,2) * t234 + t231;
t198 = Icges(3,2) * t236 + t290;
t190 = -rSges(5,3) * t248 + (rSges(5,1) * t244 - rSges(5,2) * t243) * t246;
t189 = t236 * t274 + t211;
t188 = t234 * t274 + t210;
t187 = -Icges(5,5) * t248 + (Icges(5,1) * t244 - Icges(5,4) * t243) * t246;
t186 = -Icges(5,6) * t248 + (Icges(5,4) * t244 - Icges(5,2) * t243) * t246;
t185 = -Icges(5,3) * t248 + (Icges(5,5) * t244 - Icges(5,6) * t243) * t246;
t183 = t236 * t280 + t287;
t182 = t234 * t244 - t236 * t281;
t181 = t234 * t280 - t284;
t180 = -t234 * t281 - t236 * t244;
t179 = -rSges(6,3) * t248 + (rSges(6,1) * t235 - rSges(6,2) * t233) * t246;
t165 = rSges(4,3) * t234 + t236 * t265;
t164 = -rSges(4,3) * t236 + t234 * t265;
t163 = V_base(5) * rSges(2,3) - t226 * t237 + t272;
t162 = t227 * t237 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t226 * V_base(4) - t227 * V_base(5) + V_base(3);
t152 = V_base(5) * rSges(3,3) + (-t202 - t295) * t237 + t268;
t151 = t203 * t237 + (-rSges(3,3) + t292) * V_base(4) + t273;
t150 = V_base(4) * t202 + (-t203 - t294) * V_base(5) + t267;
t145 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t283;
t144 = rSges(5,1) * t181 + rSges(5,2) * t180 + rSges(5,3) * t286;
t143 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t283;
t142 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t286;
t141 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t283;
t140 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t286;
t139 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t283;
t138 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t286;
t136 = rSges(6,1) * t171 - rSges(6,2) * t170 + rSges(6,3) * t283;
t134 = rSges(6,1) * t169 - rSges(6,2) * t168 + rSges(6,3) * t286;
t119 = t210 * t225 + (-t164 + t269) * t237 + t268;
t118 = t165 * t237 - t211 * t225 + t257;
t117 = t211 * t164 - t210 * t165 + t255;
t116 = t190 * t210 + (-t144 + t266) * t237 + t260;
t115 = t145 * t237 + (-t190 - t224) * t211 + t256;
t114 = t211 * t144 + (-t145 - t193) * t210 + t253;
t113 = -t134 * t223 + t179 * t188 + t254;
t112 = t136 * t223 - t179 * t189 + t252;
t111 = t189 * t134 - t188 * t136 + t251;
t110 = qJD(6) * t170 + t188 * t276 - t223 * t278 + t254;
t109 = qJD(6) * t168 - t189 * t276 + t223 * t277 + t252;
t108 = qJD(6) * t246 * t233 - t188 * t277 + t189 * t278 + t251;
t1 = m(1) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t150 ^ 2 + t151 ^ 2 + t152 ^ 2) / 0.2e1 + ((t168 * t306 + t169 * t304 + t286 * t305) * t223 + (t168 * t311 + t169 * t307 + t286 * t309) * t189 + (t312 * t168 + t308 * t169 + t310 * t286) * t188) * t188 / 0.2e1 + ((t170 * t306 + t171 * t304 + t283 * t305) * t223 + (t311 * t170 + t307 * t171 + t309 * t283) * t189 + (t170 * t312 + t308 * t171 + t310 * t283) * t188) * t189 / 0.2e1 + ((t139 * t286 + t141 * t180 + t143 * t181) * t211 + (t138 * t286 + t180 * t140 + t181 * t142) * t210 + (t180 * t186 + t181 * t187 + t185 * t286) * t237 + t234 * t250 - t259 * t236) * t210 / 0.2e1 + ((t139 * t283 + t182 * t141 + t183 * t143) * t211 + (t138 * t283 + t140 * t182 + t142 * t183) * t210 + (t182 * t186 + t183 * t187 + t185 * t283) * t237 + t234 * t259 + t236 * t250) * t211 / 0.2e1 + ((-t188 * t310 - t189 * t309 - t305 * t223) * t248 + ((t233 * t306 + t235 * t304) * t223 + (t311 * t233 + t307 * t235) * t189 + (t233 * t312 + t308 * t235) * t188) * t246) * t223 / 0.2e1 + ((-t198 * t234 + t200 * t236 - t247 * t216 + t219 * t249 + Icges(1,4)) * V_base(5) + (-t234 * t199 + t236 * t201 - t247 * t217 + t249 * t220 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t236 * t198 + t234 * t200 + t249 * t216 + t247 * t219 + Icges(1,2)) * V_base(5) + (t199 * t236 + t201 * t234 + t217 * t249 + t247 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t138 * t210 - t139 * t211) * t248 + ((-t141 * t243 + t143 * t244) * t211 + (-t140 * t243 + t142 * t244) * t210) * t246 + (t159 * t248 + t161 * t246) * t211 + (t158 * t248 + t160 * t246) * t210 + (Icges(2,3) + Icges(3,3) + (-t185 + t215) * t248 + (-t186 * t243 + t187 * t244 + t218) * t246) * t237) * t237 / 0.2e1 + t237 * V_base(5) * (Icges(2,5) * t247 + Icges(3,5) * t234 + Icges(2,6) * t249 + Icges(3,6) * t236) + t237 * V_base(4) * (Icges(2,5) * t249 + Icges(3,5) * t236 - Icges(2,6) * t247 - Icges(3,6) * t234) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
