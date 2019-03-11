% Calculate kinetic energy for
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:40
% EndTime: 2019-03-09 15:31:44
% DurationCPUTime: 3.53s
% Computational Cost: add. (2117->360), mult. (2971->511), div. (0->0), fcn. (3085->10), ass. (0->164)
t325 = Icges(5,1) + Icges(6,1);
t324 = -Icges(5,4) + Icges(6,5);
t323 = Icges(6,4) + Icges(5,5);
t322 = Icges(5,2) + Icges(6,3);
t321 = -Icges(6,6) + Icges(5,6);
t320 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t258 = qJ(3) + pkin(10);
t252 = sin(t258);
t253 = cos(t258);
t267 = cos(qJ(1));
t263 = sin(qJ(1));
t266 = cos(qJ(2));
t297 = t263 * t266;
t193 = t252 * t297 + t253 * t267;
t194 = -t252 * t267 + t253 * t297;
t262 = sin(qJ(2));
t300 = t262 * t263;
t319 = t193 * t322 + t194 * t324 - t300 * t321;
t296 = t266 * t267;
t195 = t252 * t296 - t263 * t253;
t196 = t263 * t252 + t253 * t296;
t299 = t262 * t267;
t318 = t195 * t322 + t196 * t324 - t299 * t321;
t317 = t324 * t193 + t194 * t325 + t323 * t300;
t316 = t324 * t195 + t196 * t325 + t323 * t299;
t315 = t321 * t266 + (t252 * t322 + t253 * t324) * t262;
t314 = -t323 * t266 + (t324 * t252 + t253 * t325) * t262;
t261 = sin(qJ(3));
t265 = cos(qJ(3));
t214 = -t261 * t297 - t265 * t267;
t301 = t261 * t267;
t215 = t265 * t297 - t301;
t313 = Icges(4,5) * t215 + Icges(4,6) * t214 - t193 * t321 + t194 * t323 - t300 * t320;
t216 = -t261 * t296 + t263 * t265;
t298 = t263 * t261;
t217 = t265 * t296 + t298;
t312 = Icges(4,5) * t217 + Icges(4,6) * t216 - t195 * t321 + t196 * t323 - t299 * t320;
t311 = t320 * t266 + (Icges(4,5) * t265 - Icges(4,6) * t261 - t252 * t321 + t253 * t323) * t262;
t306 = pkin(3) * t265;
t304 = Icges(2,4) * t263;
t303 = Icges(3,4) * t262;
t302 = Icges(3,4) * t266;
t158 = pkin(4) * t194 + qJ(5) * t193;
t279 = qJ(4) * t262 + t266 * t306;
t166 = -pkin(3) * t301 + t263 * t279;
t295 = -t158 - t166;
t159 = pkin(4) * t196 + qJ(5) * t195;
t167 = pkin(3) * t298 + t267 * t279;
t294 = -t159 - t167;
t177 = -qJ(4) * t266 + t262 * t306;
t207 = (pkin(4) * t253 + qJ(5) * t252) * t262;
t293 = -t177 - t207;
t292 = qJD(3) * t262;
t291 = qJD(4) * t262;
t290 = qJD(6) * t262;
t289 = V_base(5) * pkin(6) + V_base(1);
t246 = qJD(2) * t263 + V_base(4);
t254 = V_base(6) + qJD(1);
t213 = t267 * t292 + t246;
t286 = pkin(2) * t266 + pkin(8) * t262;
t245 = -qJD(2) * t267 + V_base(5);
t285 = rSges(3,1) * t266 - rSges(3,2) * t262;
t284 = Icges(3,1) * t266 - t303;
t283 = -Icges(3,2) * t262 + t302;
t282 = Icges(3,5) * t266 - Icges(3,6) * t262;
t212 = t263 * t292 + t245;
t243 = pkin(1) * t267 + t263 * pkin(7);
t281 = -V_base(4) * pkin(6) + t254 * t243 + V_base(2);
t242 = t263 * pkin(1) - pkin(7) * t267;
t280 = V_base(4) * t242 - t243 * V_base(5) + V_base(3);
t219 = t286 * t263;
t241 = pkin(2) * t262 - pkin(8) * t266;
t278 = t245 * t241 + (-t219 - t242) * t254 + t289;
t277 = (-Icges(3,3) * t267 + t263 * t282) * t245 + (Icges(3,3) * t263 + t267 * t282) * t246 + (Icges(3,5) * t262 + Icges(3,6) * t266) * t254;
t220 = t286 * t267;
t276 = t254 * t220 - t241 * t246 + t281;
t275 = t212 * t177 + t267 * t291 + t278;
t274 = t246 * t219 - t220 * t245 + t280;
t237 = -qJD(3) * t266 + t254;
t273 = t237 * t167 + t263 * t291 + t276;
t272 = qJD(5) * t195 + t212 * t207 + t275;
t271 = -qJD(4) * t266 + t213 * t166 + t274;
t270 = qJD(5) * t193 + t237 * t159 + t273;
t269 = qJD(5) * t262 * t252 + t213 * t158 + t271;
t202 = -Icges(3,6) * t267 + t263 * t283;
t203 = Icges(3,6) * t263 + t267 * t283;
t205 = -Icges(3,5) * t267 + t263 * t284;
t206 = Icges(3,5) * t263 + t267 * t284;
t230 = Icges(3,2) * t266 + t303;
t233 = Icges(3,1) * t262 + t302;
t268 = (-t203 * t262 + t206 * t266) * t246 + (-t202 * t262 + t205 * t266) * t245 + (-t230 * t262 + t233 * t266) * t254;
t264 = cos(qJ(6));
t260 = sin(qJ(6));
t256 = Icges(2,4) * t267;
t240 = rSges(2,1) * t267 - t263 * rSges(2,2);
t239 = t263 * rSges(2,1) + rSges(2,2) * t267;
t238 = rSges(3,1) * t262 + rSges(3,2) * t266;
t235 = Icges(2,1) * t267 - t304;
t234 = Icges(2,1) * t263 + t256;
t232 = -Icges(2,2) * t263 + t256;
t231 = Icges(2,2) * t267 + t304;
t226 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t225 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t224 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t222 = pkin(5) * t253 * t262 + pkin(9) * t266;
t221 = (-qJD(3) + qJD(6)) * t266 + t254;
t210 = t263 * rSges(3,3) + t267 * t285;
t209 = -rSges(3,3) * t267 + t263 * t285;
t208 = -rSges(4,3) * t266 + (rSges(4,1) * t265 - rSges(4,2) * t261) * t262;
t204 = -Icges(4,5) * t266 + (Icges(4,1) * t265 - Icges(4,4) * t261) * t262;
t201 = -Icges(4,6) * t266 + (Icges(4,4) * t265 - Icges(4,2) * t261) * t262;
t190 = (t252 * t260 + t253 * t264) * t262;
t189 = (t252 * t264 - t253 * t260) * t262;
t188 = -t267 * t290 + t213;
t187 = -t263 * t290 + t212;
t185 = -rSges(5,3) * t266 + (rSges(5,1) * t253 - rSges(5,2) * t252) * t262;
t184 = -rSges(6,2) * t266 + (rSges(6,1) * t253 + rSges(6,3) * t252) * t262;
t176 = V_base(5) * rSges(2,3) - t239 * t254 + t289;
t175 = t240 * t254 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t174 = t239 * V_base(4) - t240 * V_base(5) + V_base(3);
t173 = t196 * pkin(5) - pkin(9) * t299;
t172 = pkin(5) * t194 - pkin(9) * t300;
t169 = t217 * rSges(4,1) + t216 * rSges(4,2) + rSges(4,3) * t299;
t168 = rSges(4,1) * t215 + rSges(4,2) * t214 + rSges(4,3) * t300;
t165 = Icges(4,1) * t217 + Icges(4,4) * t216 + Icges(4,5) * t299;
t164 = Icges(4,1) * t215 + Icges(4,4) * t214 + Icges(4,5) * t300;
t163 = Icges(4,4) * t217 + Icges(4,2) * t216 + Icges(4,6) * t299;
t162 = Icges(4,4) * t215 + Icges(4,2) * t214 + Icges(4,6) * t300;
t157 = t195 * t260 + t196 * t264;
t156 = t195 * t264 - t196 * t260;
t155 = t193 * t260 + t194 * t264;
t154 = t193 * t264 - t194 * t260;
t153 = t196 * rSges(5,1) - t195 * rSges(5,2) + rSges(5,3) * t299;
t152 = t196 * rSges(6,1) + rSges(6,2) * t299 + t195 * rSges(6,3);
t151 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t300;
t150 = rSges(6,1) * t194 + rSges(6,2) * t300 + rSges(6,3) * t193;
t136 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t266;
t134 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t266;
t133 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t266;
t132 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t266;
t129 = t238 * t245 + (-t209 - t242) * t254 + t289;
t128 = t210 * t254 - t238 * t246 + t281;
t127 = t209 * t246 - t210 * t245 + t280;
t126 = t157 * rSges(7,1) + t156 * rSges(7,2) - rSges(7,3) * t299;
t125 = rSges(7,1) * t155 + rSges(7,2) * t154 - rSges(7,3) * t300;
t124 = Icges(7,1) * t157 + Icges(7,4) * t156 - Icges(7,5) * t299;
t123 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t300;
t122 = Icges(7,4) * t157 + Icges(7,2) * t156 - Icges(7,6) * t299;
t121 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t300;
t120 = Icges(7,5) * t157 + Icges(7,6) * t156 - Icges(7,3) * t299;
t119 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t300;
t118 = -t168 * t237 + t208 * t212 + t278;
t117 = t169 * t237 - t208 * t213 + t276;
t116 = t168 * t213 - t169 * t212 + t274;
t115 = t185 * t212 + (-t151 - t166) * t237 + t275;
t114 = t153 * t237 + (-t177 - t185) * t213 + t273;
t113 = t151 * t213 + (-t153 - t167) * t212 + t271;
t112 = t184 * t212 + (-t150 + t295) * t237 + t272;
t111 = t152 * t237 + (-t184 + t293) * t213 + t270;
t110 = t150 * t213 + (-t152 + t294) * t212 + t269;
t109 = -t125 * t221 + t136 * t187 + t212 * t222 + (-t172 + t295) * t237 + t272;
t108 = t126 * t221 - t136 * t188 + t173 * t237 + (-t222 + t293) * t213 + t270;
t107 = t125 * t188 - t126 * t187 + t172 * t213 + (-t173 + t294) * t212 + t269;
t1 = t187 * ((-t120 * t300 + t122 * t154 + t124 * t155) * t188 + (-t119 * t300 + t154 * t121 + t155 * t123) * t187 + (-t132 * t300 + t133 * t154 + t134 * t155) * t221) / 0.2e1 + t188 * ((-t120 * t299 + t156 * t122 + t157 * t124) * t188 + (-t119 * t299 + t156 * t121 + t157 * t123) * t187 + (-t132 * t299 + t156 * t133 + t157 * t134) * t221) / 0.2e1 + t221 * ((t120 * t266 + t122 * t189 + t124 * t190) * t188 + (t119 * t266 + t121 * t189 + t123 * t190) * t187 + (t266 * t132 + t189 * t133 + t190 * t134) * t221) / 0.2e1 + t246 * (t263 * t277 + t267 * t268) / 0.2e1 + t245 * (t263 * t268 - t277 * t267) / 0.2e1 + m(1) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(2) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + ((t203 * t266 + t206 * t262) * t246 + (t202 * t266 + t205 * t262) * t245 + (t266 * t230 + t262 * t233 + Icges(2,3)) * t254) * t254 / 0.2e1 + ((-t263 * t231 + t234 * t267 + Icges(1,4)) * V_base(5) + (-t263 * t232 + t267 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t231 + t263 * t234 + Icges(1,2)) * V_base(5) + (t232 * t267 + t263 * t235 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t315 * t193 + t314 * t194 + t201 * t214 + t204 * t215 + t311 * t300) * t237 + (t163 * t214 + t165 * t215 + t318 * t193 + t316 * t194 + t312 * t300) * t213 + (t214 * t162 + t215 * t164 + t319 * t193 + t317 * t194 + t313 * t300) * t212) * t212 / 0.2e1 + ((t315 * t195 + t314 * t196 + t216 * t201 + t217 * t204 + t311 * t299) * t237 + (t216 * t163 + t217 * t165 + t318 * t195 + t316 * t196 + t312 * t299) * t213 + (t216 * t162 + t217 * t164 + t319 * t195 + t317 * t196 + t313 * t299) * t212) * t213 / 0.2e1 + ((-t313 * t212 - t312 * t213 - t311 * t237) * t266 + ((-t201 * t261 + t204 * t265 + t315 * t252 + t314 * t253) * t237 + (-t163 * t261 + t165 * t265 + t318 * t252 + t316 * t253) * t213 + (-t162 * t261 + t164 * t265 + t319 * t252 + t317 * t253) * t212) * t262) * t237 / 0.2e1 + V_base(4) * t254 * (Icges(2,5) * t267 - Icges(2,6) * t263) + V_base(5) * t254 * (Icges(2,5) * t263 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
