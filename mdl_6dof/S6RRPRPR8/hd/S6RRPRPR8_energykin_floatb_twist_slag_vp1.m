% Calculate kinetic energy for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:27
% EndTime: 2019-03-09 10:49:31
% DurationCPUTime: 3.90s
% Computational Cost: add. (2084->369), mult. (2916->533), div. (0->0), fcn. (3030->10), ass. (0->165)
t325 = Icges(5,1) + Icges(6,1);
t324 = -Icges(5,4) + Icges(6,5);
t323 = Icges(6,4) + Icges(5,5);
t322 = Icges(5,2) + Icges(6,3);
t321 = -Icges(6,6) + Icges(5,6);
t320 = -Icges(5,3) - Icges(6,2);
t258 = pkin(10) + qJ(4);
t252 = sin(t258);
t253 = cos(t258);
t267 = cos(qJ(1));
t264 = sin(qJ(1));
t266 = cos(qJ(2));
t296 = t264 * t266;
t198 = t252 * t296 + t253 * t267;
t199 = -t252 * t267 + t253 * t296;
t263 = sin(qJ(2));
t299 = t263 * t264;
t319 = t198 * t322 + t199 * t324 - t299 * t321;
t295 = t266 * t267;
t200 = t252 * t295 - t264 * t253;
t201 = t264 * t252 + t253 * t295;
t298 = t263 * t267;
t318 = t200 * t322 + t201 * t324 - t298 * t321;
t317 = -t198 * t321 + t199 * t323 - t299 * t320;
t316 = -t200 * t321 + t201 * t323 - t298 * t320;
t315 = t324 * t198 + t199 * t325 + t323 * t299;
t314 = t324 * t200 + t201 * t325 + t323 * t298;
t313 = t321 * t266 + (t252 * t322 + t253 * t324) * t263;
t312 = t320 * t266 + (-t252 * t321 + t253 * t323) * t263;
t311 = -t323 * t266 + (t324 * t252 + t253 * t325) * t263;
t260 = cos(pkin(10));
t304 = pkin(3) * t260;
t303 = Icges(2,4) * t264;
t302 = Icges(3,4) * t263;
t301 = Icges(3,4) * t266;
t259 = sin(pkin(10));
t300 = t259 * t267;
t297 = t264 * t259;
t284 = pkin(2) * t266 + qJ(3) * t263;
t218 = t284 * t264;
t242 = t264 * pkin(1) - pkin(7) * t267;
t293 = -t218 - t242;
t292 = qJD(3) * t263;
t291 = qJD(4) * t263;
t290 = qJD(6) * t263;
t289 = V_base(5) * pkin(6) + V_base(1);
t246 = qJD(2) * t264 + V_base(4);
t254 = V_base(6) + qJD(1);
t217 = t267 * t291 + t246;
t238 = pkin(2) * t263 - qJ(3) * t266;
t245 = -qJD(2) * t267 + V_base(5);
t286 = t245 * t238 + t267 * t292 + t289;
t285 = rSges(3,1) * t266 - rSges(3,2) * t263;
t283 = Icges(3,1) * t266 - t302;
t282 = -Icges(3,2) * t263 + t301;
t281 = Icges(3,5) * t266 - Icges(3,6) * t263;
t216 = t264 * t291 + t245;
t243 = pkin(1) * t267 + t264 * pkin(7);
t280 = -V_base(4) * pkin(6) + t254 * t243 + V_base(2);
t279 = V_base(4) * t242 - t243 * V_base(5) + V_base(3);
t278 = (-Icges(3,3) * t267 + t264 * t281) * t245 + (Icges(3,3) * t264 + t267 * t281) * t246 + (Icges(3,5) * t263 + Icges(3,6) * t266) * t254;
t277 = pkin(8) * t263 + t266 * t304;
t219 = t284 * t267;
t276 = t254 * t219 + t264 * t292 + t280;
t275 = -qJD(3) * t266 + t246 * t218 + t279;
t168 = -pkin(3) * t300 + t264 * t277;
t177 = -pkin(8) * t266 + t263 * t304;
t274 = t245 * t177 + (-t168 + t293) * t254 + t286;
t208 = (pkin(4) * t253 + qJ(5) * t252) * t263;
t273 = qJD(5) * t200 + t216 * t208 + t274;
t169 = pkin(3) * t297 + t267 * t277;
t272 = t254 * t169 + (-t177 - t238) * t246 + t276;
t271 = t246 * t168 + (-t169 - t219) * t245 + t275;
t159 = pkin(4) * t201 + qJ(5) * t200;
t237 = -qJD(4) * t266 + t254;
t270 = qJD(5) * t198 + t237 * t159 + t272;
t158 = pkin(4) * t199 + qJ(5) * t198;
t269 = qJD(5) * t263 * t252 + t217 * t158 + t271;
t204 = -Icges(3,6) * t267 + t264 * t282;
t205 = Icges(3,6) * t264 + t267 * t282;
t206 = -Icges(3,5) * t267 + t264 * t283;
t207 = Icges(3,5) * t264 + t267 * t283;
t230 = Icges(3,2) * t266 + t302;
t233 = Icges(3,1) * t263 + t301;
t268 = (-t205 * t263 + t207 * t266) * t246 + (-t204 * t263 + t206 * t266) * t245 + (-t230 * t263 + t233 * t266) * t254;
t265 = cos(qJ(6));
t262 = sin(qJ(6));
t256 = Icges(2,4) * t267;
t241 = rSges(2,1) * t267 - t264 * rSges(2,2);
t240 = t264 * rSges(2,1) + rSges(2,2) * t267;
t239 = rSges(3,1) * t263 + rSges(3,2) * t266;
t235 = Icges(2,1) * t267 - t303;
t234 = Icges(2,1) * t264 + t256;
t232 = -Icges(2,2) * t264 + t256;
t231 = Icges(2,2) * t267 + t303;
t226 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t225 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t224 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t222 = pkin(5) * t253 * t263 + pkin(9) * t266;
t221 = (-qJD(4) + qJD(6)) * t266 + t254;
t215 = t260 * t295 + t297;
t214 = -t259 * t295 + t264 * t260;
t213 = t260 * t296 - t300;
t212 = -t259 * t296 - t260 * t267;
t210 = t264 * rSges(3,3) + t267 * t285;
t209 = -rSges(3,3) * t267 + t264 * t285;
t197 = -rSges(4,3) * t266 + (rSges(4,1) * t260 - rSges(4,2) * t259) * t263;
t195 = -Icges(4,5) * t266 + (Icges(4,1) * t260 - Icges(4,4) * t259) * t263;
t194 = -Icges(4,6) * t266 + (Icges(4,4) * t260 - Icges(4,2) * t259) * t263;
t193 = -Icges(4,3) * t266 + (Icges(4,5) * t260 - Icges(4,6) * t259) * t263;
t190 = (t252 * t262 + t253 * t265) * t263;
t189 = (t252 * t265 - t253 * t262) * t263;
t188 = -t267 * t290 + t217;
t187 = -t264 * t290 + t216;
t185 = -rSges(5,3) * t266 + (rSges(5,1) * t253 - rSges(5,2) * t252) * t263;
t184 = -rSges(6,2) * t266 + (rSges(6,1) * t253 + rSges(6,3) * t252) * t263;
t176 = V_base(5) * rSges(2,3) - t240 * t254 + t289;
t175 = t241 * t254 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t174 = t240 * V_base(4) - t241 * V_base(5) + V_base(3);
t173 = t201 * pkin(5) - pkin(9) * t298;
t172 = pkin(5) * t199 - pkin(9) * t299;
t167 = t215 * rSges(4,1) + t214 * rSges(4,2) + rSges(4,3) * t298;
t166 = rSges(4,1) * t213 + rSges(4,2) * t212 + rSges(4,3) * t299;
t165 = Icges(4,1) * t215 + Icges(4,4) * t214 + Icges(4,5) * t298;
t164 = Icges(4,1) * t213 + Icges(4,4) * t212 + Icges(4,5) * t299;
t163 = Icges(4,4) * t215 + Icges(4,2) * t214 + Icges(4,6) * t298;
t162 = Icges(4,4) * t213 + Icges(4,2) * t212 + Icges(4,6) * t299;
t161 = Icges(4,5) * t215 + Icges(4,6) * t214 + Icges(4,3) * t298;
t160 = Icges(4,5) * t213 + Icges(4,6) * t212 + Icges(4,3) * t299;
t157 = t200 * t262 + t201 * t265;
t156 = t200 * t265 - t201 * t262;
t155 = t198 * t262 + t199 * t265;
t154 = t198 * t265 - t199 * t262;
t152 = t201 * rSges(5,1) - t200 * rSges(5,2) + rSges(5,3) * t298;
t151 = t201 * rSges(6,1) + rSges(6,2) * t298 + t200 * rSges(6,3);
t150 = rSges(5,1) * t199 - rSges(5,2) * t198 + rSges(5,3) * t299;
t149 = rSges(6,1) * t199 + rSges(6,2) * t299 + rSges(6,3) * t198;
t135 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t266;
t133 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t266;
t132 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t266;
t131 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t266;
t129 = t239 * t245 + (-t209 - t242) * t254 + t289;
t128 = t210 * t254 - t239 * t246 + t280;
t127 = t209 * t246 - t210 * t245 + t279;
t126 = t157 * rSges(7,1) + t156 * rSges(7,2) - rSges(7,3) * t298;
t125 = rSges(7,1) * t155 + rSges(7,2) * t154 - rSges(7,3) * t299;
t124 = Icges(7,1) * t157 + Icges(7,4) * t156 - Icges(7,5) * t298;
t123 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t299;
t122 = Icges(7,4) * t157 + Icges(7,2) * t156 - Icges(7,6) * t298;
t121 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t299;
t120 = Icges(7,5) * t157 + Icges(7,6) * t156 - Icges(7,3) * t298;
t119 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t299;
t118 = t197 * t245 + (-t166 + t293) * t254 + t286;
t117 = t167 * t254 + (-t197 - t238) * t246 + t276;
t116 = t166 * t246 + (-t167 - t219) * t245 + t275;
t115 = -t150 * t237 + t185 * t216 + t274;
t114 = t152 * t237 - t185 * t217 + t272;
t113 = t150 * t217 - t152 * t216 + t271;
t112 = t184 * t216 + (-t149 - t158) * t237 + t273;
t111 = t151 * t237 + (-t184 - t208) * t217 + t270;
t110 = t149 * t217 + (-t151 - t159) * t216 + t269;
t109 = -t125 * t221 + t135 * t187 + t216 * t222 + (-t158 - t172) * t237 + t273;
t108 = t126 * t221 - t135 * t188 + t173 * t237 + (-t208 - t222) * t217 + t270;
t107 = t125 * t188 - t126 * t187 + t172 * t217 + (-t159 - t173) * t216 + t269;
t1 = t187 * ((-t120 * t299 + t122 * t154 + t124 * t155) * t188 + (-t119 * t299 + t154 * t121 + t155 * t123) * t187 + (-t131 * t299 + t132 * t154 + t133 * t155) * t221) / 0.2e1 + t188 * ((-t120 * t298 + t156 * t122 + t157 * t124) * t188 + (-t119 * t298 + t156 * t121 + t157 * t123) * t187 + (-t131 * t298 + t156 * t132 + t157 * t133) * t221) / 0.2e1 + m(2) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t221 * ((t120 * t266 + t122 * t189 + t124 * t190) * t188 + (t119 * t266 + t121 * t189 + t123 * t190) * t187 + (t266 * t131 + t189 * t132 + t190 * t133) * t221) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(1) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + ((t198 * t313 + t199 * t311 + t299 * t312) * t237 + (t198 * t318 + t199 * t314 + t299 * t316) * t217 + (t319 * t198 + t315 * t199 + t317 * t299) * t216) * t216 / 0.2e1 + ((t200 * t313 + t201 * t311 + t298 * t312) * t237 + (t318 * t200 + t314 * t201 + t316 * t298) * t217 + (t200 * t319 + t315 * t201 + t317 * t298) * t216) * t217 / 0.2e1 + ((-t216 * t317 - t217 * t316 - t312 * t237) * t266 + ((t313 * t252 + t311 * t253) * t237 + (t252 * t318 + t253 * t314) * t217 + (t252 * t319 + t315 * t253) * t216) * t263) * t237 / 0.2e1 + ((t161 * t299 + t163 * t212 + t165 * t213) * t246 + (t160 * t299 + t212 * t162 + t213 * t164) * t245 + (t193 * t299 + t194 * t212 + t195 * t213) * t254 + t264 * t268 - t267 * t278) * t245 / 0.2e1 + ((t161 * t298 + t214 * t163 + t215 * t165) * t246 + (t160 * t298 + t214 * t162 + t215 * t164) * t245 + (t193 * t298 + t214 * t194 + t215 * t195) * t254 + t264 * t278 + t267 * t268) * t246 / 0.2e1 + ((-t264 * t231 + t234 * t267 + Icges(1,4)) * V_base(5) + (-t264 * t232 + t267 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t231 + t264 * t234 + Icges(1,2)) * V_base(5) + (t232 * t267 + t264 * t235 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t160 * t245 - t161 * t246) * t266 + ((-t163 * t259 + t165 * t260) * t246 + (-t162 * t259 + t164 * t260) * t245) * t263 + (t205 * t266 + t207 * t263) * t246 + (t204 * t266 + t206 * t263) * t245 + (Icges(2,3) + (-t193 + t230) * t266 + (-t194 * t259 + t195 * t260 + t233) * t263) * t254) * t254 / 0.2e1 + t254 * V_base(4) * (Icges(2,5) * t267 - Icges(2,6) * t264) + t254 * V_base(5) * (Icges(2,5) * t264 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
