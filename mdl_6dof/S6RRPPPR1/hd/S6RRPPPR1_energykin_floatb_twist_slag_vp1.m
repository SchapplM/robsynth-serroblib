% Calculate kinetic energy for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:05:59
% EndTime: 2019-03-09 08:06:02
% DurationCPUTime: 2.92s
% Computational Cost: add. (1939->342), mult. (2511->463), div. (0->0), fcn. (2569->10), ass. (0->165)
t328 = Icges(5,1) + Icges(6,1);
t327 = -Icges(5,4) + Icges(6,5);
t326 = Icges(6,4) + Icges(5,5);
t325 = Icges(5,2) + Icges(6,3);
t324 = -Icges(6,6) + Icges(5,6);
t323 = Icges(3,3) + Icges(4,3);
t322 = -Icges(5,3) - Icges(6,2);
t247 = qJ(2) + pkin(9);
t240 = sin(t247);
t241 = cos(t247);
t252 = sin(qJ(2));
t255 = cos(qJ(2));
t321 = Icges(3,5) * t255 + Icges(4,5) * t241 - Icges(3,6) * t252 - Icges(4,6) * t240;
t248 = sin(pkin(10));
t253 = sin(qJ(1));
t294 = t253 * t248;
t249 = cos(pkin(10));
t256 = cos(qJ(1));
t295 = t249 * t256;
t200 = t241 * t294 + t295;
t293 = t253 * t249;
t296 = t248 * t256;
t201 = t241 * t293 - t296;
t298 = t240 * t253;
t320 = t325 * t200 + t327 * t201 - t324 * t298;
t202 = t241 * t296 - t293;
t203 = t241 * t295 + t294;
t297 = t240 * t256;
t319 = t325 * t202 + t327 * t203 - t324 * t297;
t318 = -t324 * t200 + t326 * t201 - t322 * t298;
t317 = -t324 * t202 + t326 * t203 - t322 * t297;
t316 = t327 * t200 + t328 * t201 + t326 * t298;
t315 = t327 * t202 + t328 * t203 + t326 * t297;
t314 = t324 * t241 + (t325 * t248 + t327 * t249) * t240;
t313 = t322 * t241 + (-t324 * t248 + t326 * t249) * t240;
t312 = -t326 * t241 + (t327 * t248 + t328 * t249) * t240;
t299 = Icges(4,4) * t241;
t273 = -Icges(4,2) * t240 + t299;
t178 = -Icges(4,6) * t256 + t253 * t273;
t179 = Icges(4,6) * t253 + t256 * t273;
t300 = Icges(4,4) * t240;
t275 = Icges(4,1) * t241 - t300;
t180 = -Icges(4,5) * t256 + t253 * t275;
t181 = Icges(4,5) * t253 + t256 * t275;
t301 = Icges(3,4) * t255;
t274 = -Icges(3,2) * t252 + t301;
t192 = -Icges(3,6) * t256 + t253 * t274;
t193 = Icges(3,6) * t253 + t256 * t274;
t302 = Icges(3,4) * t252;
t276 = Icges(3,1) * t255 - t302;
t194 = -Icges(3,5) * t256 + t253 * t276;
t195 = Icges(3,5) * t253 + t256 * t276;
t209 = Icges(4,2) * t241 + t300;
t210 = Icges(4,1) * t240 + t299;
t223 = Icges(3,2) * t255 + t302;
t226 = Icges(3,1) * t252 + t301;
t236 = -qJD(2) * t256 + V_base(5);
t237 = qJD(2) * t253 + V_base(4);
t242 = V_base(6) + qJD(1);
t311 = (-t209 * t240 + t210 * t241 - t223 * t252 + t226 * t255) * t242 + (-t179 * t240 + t181 * t241 - t193 * t252 + t195 * t255) * t237 + (-t178 * t240 + t180 * t241 - t192 * t252 + t194 * t255) * t236;
t310 = (Icges(3,5) * t252 + Icges(4,5) * t240 + Icges(3,6) * t255 + Icges(4,6) * t241) * t242 + (t323 * t253 + t321 * t256) * t237 + (t321 * t253 - t323 * t256) * t236;
t306 = pkin(2) * t252;
t305 = pkin(2) * t255;
t303 = Icges(2,4) * t253;
t174 = -qJ(3) * t256 + t253 * t305;
t234 = t253 * pkin(1) - pkin(7) * t256;
t292 = -t174 - t234;
t175 = qJ(3) * t253 + t256 * t305;
t277 = pkin(3) * t241 + qJ(4) * t240;
t197 = t277 * t256;
t291 = -t175 - t197;
t290 = qJD(4) * t240;
t289 = qJD(6) * t240;
t288 = V_base(5) * pkin(6) + V_base(1);
t155 = pkin(4) * t203 + qJ(5) * t202;
t285 = -t155 + t291;
t196 = t277 * t253;
t284 = -t196 + t292;
t211 = pkin(3) * t240 - qJ(4) * t241;
t283 = -t211 - t306;
t154 = pkin(4) * t201 + qJ(5) * t200;
t282 = -t154 + t284;
t189 = (pkin(4) * t249 + qJ(5) * t248) * t240;
t281 = -t189 + t283;
t280 = qJD(3) * t253 + t236 * t306 + t288;
t279 = rSges(3,1) * t255 - rSges(3,2) * t252;
t278 = rSges(4,1) * t241 - rSges(4,2) * t240;
t235 = pkin(1) * t256 + t253 * pkin(7);
t270 = -V_base(4) * pkin(6) + t242 * t235 + V_base(2);
t269 = V_base(4) * t234 - t235 * V_base(5) + V_base(3);
t268 = t236 * t211 + t256 * t290 + t280;
t267 = t237 * t174 + t269;
t264 = qJD(5) * t202 + t236 * t189 + t268;
t263 = -qJD(3) * t256 + t242 * t175 + t270;
t262 = -qJD(4) * t241 + t237 * t196 + t267;
t261 = t242 * t197 + t253 * t290 + t263;
t260 = qJD(5) * t240 * t248 + t237 * t154 + t262;
t259 = qJD(5) * t200 + t242 * t155 + t261;
t254 = cos(qJ(6));
t251 = sin(qJ(6));
t245 = Icges(2,4) * t256;
t233 = rSges(2,1) * t256 - t253 * rSges(2,2);
t232 = t253 * rSges(2,1) + rSges(2,2) * t256;
t231 = rSges(3,1) * t252 + rSges(3,2) * t255;
t228 = Icges(2,1) * t256 - t303;
t227 = Icges(2,1) * t253 + t245;
t225 = -Icges(2,2) * t253 + t245;
t224 = Icges(2,2) * t256 + t303;
t218 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t217 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t216 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t215 = qJD(6) * t241 + t242;
t212 = rSges(4,1) * t240 + rSges(4,2) * t241;
t206 = pkin(5) * t240 * t249 + pkin(8) * t241;
t205 = -t256 * t289 + t237;
t204 = -t253 * t289 + t236;
t199 = t253 * rSges(3,3) + t256 * t279;
t198 = -rSges(3,3) * t256 + t253 * t279;
t185 = (t248 * t251 + t249 * t254) * t240;
t184 = (t248 * t254 - t249 * t251) * t240;
t183 = t253 * rSges(4,3) + t256 * t278;
t182 = -rSges(4,3) * t256 + t253 * t278;
t172 = V_base(5) * rSges(2,3) - t232 * t242 + t288;
t171 = t233 * t242 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t170 = -rSges(5,3) * t241 + (rSges(5,1) * t249 - rSges(5,2) * t248) * t240;
t169 = -rSges(6,2) * t241 + (rSges(6,1) * t249 + rSges(6,3) * t248) * t240;
t161 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t159 = t203 * pkin(5) - pkin(8) * t297;
t158 = pkin(5) * t201 - pkin(8) * t298;
t153 = t202 * t251 + t203 * t254;
t152 = t202 * t254 - t203 * t251;
t151 = t200 * t251 + t201 * t254;
t150 = t200 * t254 - t201 * t251;
t147 = t203 * rSges(5,1) - t202 * rSges(5,2) + rSges(5,3) * t297;
t146 = t203 * rSges(6,1) + rSges(6,2) * t297 + t202 * rSges(6,3);
t145 = rSges(5,1) * t201 - rSges(5,2) * t200 + rSges(5,3) * t298;
t144 = rSges(6,1) * t201 + rSges(6,2) * t298 + rSges(6,3) * t200;
t131 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t241;
t130 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t241;
t129 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t241;
t128 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t241;
t127 = t231 * t236 + (-t198 - t234) * t242 + t288;
t126 = t199 * t242 - t231 * t237 + t270;
t125 = t198 * t237 - t199 * t236 + t269;
t124 = t153 * rSges(7,1) + t152 * rSges(7,2) - rSges(7,3) * t297;
t123 = rSges(7,1) * t151 + rSges(7,2) * t150 - rSges(7,3) * t298;
t122 = Icges(7,1) * t153 + Icges(7,4) * t152 - Icges(7,5) * t297;
t121 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t298;
t120 = Icges(7,4) * t153 + Icges(7,2) * t152 - Icges(7,6) * t297;
t119 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t298;
t118 = Icges(7,5) * t153 + Icges(7,6) * t152 - Icges(7,3) * t297;
t117 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t298;
t116 = t212 * t236 + (-t182 + t292) * t242 + t280;
t115 = t242 * t183 + (-t212 - t306) * t237 + t263;
t114 = t182 * t237 + (-t175 - t183) * t236 + t267;
t113 = t170 * t236 + (-t145 + t284) * t242 + t268;
t112 = t242 * t147 + (-t170 + t283) * t237 + t261;
t111 = t145 * t237 + (-t147 + t291) * t236 + t262;
t110 = t169 * t236 + (-t144 + t282) * t242 + t264;
t109 = t242 * t146 + (-t169 + t281) * t237 + t259;
t108 = t144 * t237 + (-t146 + t285) * t236 + t260;
t107 = -t123 * t215 + t131 * t204 + t206 * t236 + (-t158 + t282) * t242 + t264;
t106 = t215 * t124 - t205 * t131 + t242 * t159 + (-t206 + t281) * t237 + t259;
t105 = t123 * t205 - t124 * t204 + t158 * t237 + (-t159 + t285) * t236 + t260;
t1 = t215 * ((t118 * t241 + t120 * t184 + t122 * t185) * t205 + (t117 * t241 + t119 * t184 + t121 * t185) * t204 + (t241 * t128 + t184 * t129 + t185 * t130) * t215) / 0.2e1 + m(1) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + m(2) * (t161 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(3) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(6) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + t205 * ((-t118 * t297 + t152 * t120 + t153 * t122) * t205 + (-t117 * t297 + t152 * t119 + t153 * t121) * t204 + (-t128 * t297 + t152 * t129 + t153 * t130) * t215) / 0.2e1 + t204 * ((-t118 * t298 + t120 * t150 + t122 * t151) * t205 + (-t117 * t298 + t150 * t119 + t151 * t121) * t204 + (-t128 * t298 + t129 * t150 + t130 * t151) * t215) / 0.2e1 + ((-t253 * t224 + t227 * t256 + Icges(1,4)) * V_base(5) + (-t253 * t225 + t228 * t256 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t224 * t256 + t253 * t227 + Icges(1,2)) * V_base(5) + (t225 * t256 + t253 * t228 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (-t310 * t256 + t311 * t253 + (t314 * t200 + t312 * t201 + t313 * t298) * t242 + (t319 * t200 + t315 * t201 + t317 * t298) * t237 + (t320 * t200 + t316 * t201 + t318 * t298) * t236) * t236 / 0.2e1 + (t311 * t256 + t310 * t253 + (t314 * t202 + t312 * t203 + t313 * t297) * t242 + (t319 * t202 + t315 * t203 + t317 * t297) * t237 + (t320 * t202 + t316 * t203 + t318 * t297) * t236) * t237 / 0.2e1 + ((t193 * t255 + t195 * t252 + (t179 - t317) * t241 + (t319 * t248 + t315 * t249 + t181) * t240) * t237 + (t192 * t255 + t194 * t252 + (t178 - t318) * t241 + (t320 * t248 + t316 * t249 + t180) * t240) * t236 + (t223 * t255 + t226 * t252 + Icges(2,3) + (t209 - t313) * t241 + (t314 * t248 + t312 * t249 + t210) * t240) * t242) * t242 / 0.2e1 + t242 * V_base(4) * (Icges(2,5) * t256 - Icges(2,6) * t253) + t242 * V_base(5) * (Icges(2,5) * t253 + Icges(2,6) * t256) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
