% Calculate kinetic energy for
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:45
% EndTime: 2019-03-09 22:38:49
% DurationCPUTime: 3.99s
% Computational Cost: add. (2180->367), mult. (3076->544), div. (0->0), fcn. (3190->10), ass. (0->168)
t326 = Icges(5,1) + Icges(6,1);
t325 = -Icges(5,4) + Icges(6,5);
t324 = Icges(6,4) + Icges(5,5);
t323 = Icges(5,2) + Icges(6,3);
t322 = -Icges(6,6) + Icges(5,6);
t321 = -Icges(5,3) - Icges(6,2);
t261 = qJ(3) + qJ(4);
t257 = sin(t261);
t258 = cos(t261);
t269 = cos(qJ(1));
t265 = sin(qJ(1));
t268 = cos(qJ(2));
t298 = t265 * t268;
t208 = t257 * t298 + t258 * t269;
t209 = -t257 * t269 + t258 * t298;
t264 = sin(qJ(2));
t300 = t264 * t265;
t320 = t323 * t208 + t325 * t209 - t322 * t300;
t297 = t268 * t269;
t210 = t257 * t297 - t265 * t258;
t211 = t257 * t265 + t258 * t297;
t299 = t264 * t269;
t319 = t323 * t210 + t325 * t211 - t322 * t299;
t318 = -t322 * t208 + t324 * t209 - t321 * t300;
t317 = -t322 * t210 + t324 * t211 - t321 * t299;
t316 = t325 * t208 + t326 * t209 + t324 * t300;
t315 = t325 * t210 + t326 * t211 + t324 * t299;
t314 = t322 * t268 + (t323 * t257 + t325 * t258) * t264;
t313 = t321 * t268 + (-t322 * t257 + t324 * t258) * t264;
t312 = -t324 * t268 + (t325 * t257 + t326 * t258) * t264;
t267 = cos(qJ(3));
t307 = pkin(3) * t267;
t305 = Icges(2,4) * t265;
t304 = Icges(3,4) * t264;
t303 = Icges(3,4) * t268;
t263 = sin(qJ(3));
t302 = t263 * t265;
t301 = t263 * t269;
t296 = qJD(3) * t264;
t295 = qJD(4) * t264;
t294 = qJD(6) * t264;
t293 = -qJD(3) - qJD(4);
t292 = V_base(5) * pkin(6) + V_base(1);
t248 = qJD(2) * t265 + V_base(4);
t255 = V_base(6) + qJD(1);
t216 = t269 * t296 + t248;
t289 = pkin(2) * t268 + pkin(8) * t264;
t247 = -qJD(2) * t269 + V_base(5);
t288 = rSges(3,1) * t268 - rSges(3,2) * t264;
t190 = t269 * t295 + t216;
t287 = Icges(3,1) * t268 - t304;
t286 = -Icges(3,2) * t264 + t303;
t285 = Icges(3,5) * t268 - Icges(3,6) * t264;
t215 = t265 * t296 + t247;
t246 = pkin(1) * t269 + pkin(7) * t265;
t284 = -V_base(4) * pkin(6) + t255 * t246 + V_base(2);
t245 = pkin(1) * t265 - pkin(7) * t269;
t283 = V_base(4) * t245 - t246 * V_base(5) + V_base(3);
t189 = t265 * t295 + t215;
t282 = pkin(9) * t264 + t268 * t307;
t222 = t289 * t265;
t244 = t264 * pkin(2) - t268 * pkin(8);
t281 = t247 * t244 + (-t222 - t245) * t255 + t292;
t280 = (-Icges(3,3) * t269 + t265 * t285) * t247 + (Icges(3,3) * t265 + t269 * t285) * t248 + (Icges(3,5) * t264 + Icges(3,6) * t268) * t255;
t223 = t289 * t269;
t279 = t255 * t223 - t244 * t248 + t284;
t278 = t248 * t222 - t223 * t247 + t283;
t169 = -pkin(3) * t301 + t265 * t282;
t179 = -pkin(9) * t268 + t264 * t307;
t239 = -qJD(3) * t268 + t255;
t277 = -t169 * t239 + t215 * t179 + t281;
t170 = pkin(3) * t302 + t269 * t282;
t276 = t239 * t170 - t179 * t216 + t279;
t213 = (pkin(4) * t258 + qJ(5) * t257) * t264;
t275 = qJD(5) * t210 + t189 * t213 + t277;
t274 = t216 * t169 - t170 * t215 + t278;
t159 = pkin(4) * t211 + qJ(5) * t210;
t224 = t268 * t293 + t255;
t273 = qJD(5) * t208 + t224 * t159 + t276;
t158 = pkin(4) * t209 + qJ(5) * t208;
t272 = qJD(5) * t264 * t257 + t190 * t158 + t274;
t200 = -Icges(3,6) * t269 + t265 * t286;
t201 = Icges(3,6) * t265 + t269 * t286;
t203 = -Icges(3,5) * t269 + t265 * t287;
t204 = Icges(3,5) * t265 + t269 * t287;
t233 = Icges(3,2) * t268 + t304;
t236 = Icges(3,1) * t264 + t303;
t271 = (-t201 * t264 + t204 * t268) * t248 + (-t200 * t264 + t203 * t268) * t247 + (-t233 * t264 + t236 * t268) * t255;
t266 = cos(qJ(6));
t262 = sin(qJ(6));
t259 = Icges(2,4) * t269;
t242 = rSges(2,1) * t269 - rSges(2,2) * t265;
t241 = rSges(2,1) * t265 + rSges(2,2) * t269;
t240 = rSges(3,1) * t264 + rSges(3,2) * t268;
t238 = Icges(2,1) * t269 - t305;
t237 = Icges(2,1) * t265 + t259;
t235 = -Icges(2,2) * t265 + t259;
t234 = Icges(2,2) * t269 + t305;
t229 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t228 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t227 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t225 = pkin(5) * t258 * t264 + pkin(10) * t268;
t220 = t267 * t297 + t302;
t219 = -t263 * t297 + t265 * t267;
t218 = t267 * t298 - t301;
t217 = -t263 * t298 - t267 * t269;
t212 = (qJD(6) + t293) * t268 + t255;
t207 = rSges(3,3) * t265 + t269 * t288;
t206 = -rSges(3,3) * t269 + t265 * t288;
t205 = -rSges(4,3) * t268 + (rSges(4,1) * t267 - rSges(4,2) * t263) * t264;
t202 = -Icges(4,5) * t268 + (Icges(4,1) * t267 - Icges(4,4) * t263) * t264;
t199 = -Icges(4,6) * t268 + (Icges(4,4) * t267 - Icges(4,2) * t263) * t264;
t196 = -Icges(4,3) * t268 + (Icges(4,5) * t267 - Icges(4,6) * t263) * t264;
t192 = (t257 * t262 + t258 * t266) * t264;
t191 = (t257 * t266 - t258 * t262) * t264;
t188 = -rSges(5,3) * t268 + (rSges(5,1) * t258 - rSges(5,2) * t257) * t264;
t187 = -rSges(6,2) * t268 + (rSges(6,1) * t258 + rSges(6,3) * t257) * t264;
t178 = V_base(5) * rSges(2,3) - t241 * t255 + t292;
t177 = t242 * t255 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t176 = pkin(5) * t211 - pkin(10) * t299;
t175 = pkin(5) * t209 - pkin(10) * t300;
t174 = t241 * V_base(4) - t242 * V_base(5) + V_base(3);
t173 = -t269 * t294 + t190;
t172 = -t265 * t294 + t189;
t168 = rSges(4,1) * t220 + rSges(4,2) * t219 + rSges(4,3) * t299;
t167 = rSges(4,1) * t218 + rSges(4,2) * t217 + rSges(4,3) * t300;
t165 = Icges(4,1) * t220 + Icges(4,4) * t219 + Icges(4,5) * t299;
t164 = Icges(4,1) * t218 + Icges(4,4) * t217 + Icges(4,5) * t300;
t163 = Icges(4,4) * t220 + Icges(4,2) * t219 + Icges(4,6) * t299;
t162 = Icges(4,4) * t218 + Icges(4,2) * t217 + Icges(4,6) * t300;
t161 = Icges(4,5) * t220 + Icges(4,6) * t219 + Icges(4,3) * t299;
t160 = Icges(4,5) * t218 + Icges(4,6) * t217 + Icges(4,3) * t300;
t157 = t210 * t262 + t211 * t266;
t156 = t210 * t266 - t211 * t262;
t155 = t208 * t262 + t209 * t266;
t154 = t208 * t266 - t209 * t262;
t153 = rSges(5,1) * t211 - rSges(5,2) * t210 + rSges(5,3) * t299;
t152 = rSges(6,1) * t211 + rSges(6,2) * t299 + rSges(6,3) * t210;
t151 = rSges(5,1) * t209 - rSges(5,2) * t208 + rSges(5,3) * t300;
t150 = rSges(6,1) * t209 + rSges(6,2) * t300 + rSges(6,3) * t208;
t136 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t268;
t135 = Icges(7,1) * t192 + Icges(7,4) * t191 + Icges(7,5) * t268;
t134 = Icges(7,4) * t192 + Icges(7,2) * t191 + Icges(7,6) * t268;
t133 = Icges(7,5) * t192 + Icges(7,6) * t191 + Icges(7,3) * t268;
t130 = t240 * t247 + (-t206 - t245) * t255 + t292;
t129 = t207 * t255 - t240 * t248 + t284;
t127 = t206 * t248 - t207 * t247 + t283;
t126 = rSges(7,1) * t157 + rSges(7,2) * t156 - rSges(7,3) * t299;
t125 = rSges(7,1) * t155 + rSges(7,2) * t154 - rSges(7,3) * t300;
t124 = Icges(7,1) * t157 + Icges(7,4) * t156 - Icges(7,5) * t299;
t123 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t300;
t122 = Icges(7,4) * t157 + Icges(7,2) * t156 - Icges(7,6) * t299;
t121 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t300;
t120 = Icges(7,5) * t157 + Icges(7,6) * t156 - Icges(7,3) * t299;
t119 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t300;
t118 = -t167 * t239 + t205 * t215 + t281;
t117 = t168 * t239 - t205 * t216 + t279;
t116 = t167 * t216 - t168 * t215 + t278;
t115 = -t151 * t224 + t188 * t189 + t277;
t114 = t153 * t224 - t188 * t190 + t276;
t113 = t151 * t190 - t153 * t189 + t274;
t112 = t187 * t189 + (-t150 - t158) * t224 + t275;
t111 = t152 * t224 + (-t187 - t213) * t190 + t273;
t110 = t150 * t190 + (-t152 - t159) * t189 + t272;
t109 = -t125 * t212 + t136 * t172 + t189 * t225 + (-t158 - t175) * t224 + t275;
t108 = t126 * t212 - t136 * t173 + t176 * t224 + (-t213 - t225) * t190 + t273;
t107 = t125 * t173 - t126 * t172 + t175 * t190 + (-t159 - t176) * t189 + t272;
t1 = t239 * ((-t160 * t215 - t161 * t216 - t196 * t239) * t268 + ((-t163 * t263 + t165 * t267) * t216 + (-t162 * t263 + t164 * t267) * t215 + (-t199 * t263 + t202 * t267) * t239) * t264) / 0.2e1 + t215 * ((t161 * t300 + t163 * t217 + t165 * t218) * t216 + (t160 * t300 + t217 * t162 + t218 * t164) * t215 + (t196 * t300 + t199 * t217 + t202 * t218) * t239) / 0.2e1 + t173 * ((-t120 * t299 + t156 * t122 + t157 * t124) * t173 + (-t119 * t299 + t121 * t156 + t123 * t157) * t172 + (-t133 * t299 + t134 * t156 + t135 * t157) * t212) / 0.2e1 + t216 * ((t161 * t299 + t219 * t163 + t220 * t165) * t216 + (t160 * t299 + t162 * t219 + t164 * t220) * t215 + (t196 * t299 + t199 * t219 + t202 * t220) * t239) / 0.2e1 + t172 * ((-t120 * t300 + t122 * t154 + t124 * t155) * t173 + (-t119 * t300 + t154 * t121 + t155 * t123) * t172 + (-t133 * t300 + t134 * t154 + t135 * t155) * t212) / 0.2e1 + t212 * ((t120 * t268 + t122 * t191 + t124 * t192) * t173 + (t119 * t268 + t121 * t191 + t123 * t192) * t172 + (t268 * t133 + t191 * t134 + t192 * t135) * t212) / 0.2e1 + m(1) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + m(2) * (t174 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t248 * (t265 * t280 + t269 * t271) / 0.2e1 + t247 * (t265 * t271 - t280 * t269) / 0.2e1 + ((t208 * t314 + t209 * t312 + t300 * t313) * t224 + (t208 * t319 + t209 * t315 + t300 * t317) * t190 + (t320 * t208 + t316 * t209 + t318 * t300) * t189) * t189 / 0.2e1 + ((t210 * t314 + t211 * t312 + t299 * t313) * t224 + (t319 * t210 + t315 * t211 + t317 * t299) * t190 + (t210 * t320 + t316 * t211 + t318 * t299) * t189) * t190 / 0.2e1 + ((-t189 * t318 - t190 * t317 - t224 * t313) * t268 + ((t257 * t314 + t258 * t312) * t224 + (t257 * t319 + t258 * t315) * t190 + (t257 * t320 + t316 * t258) * t189) * t264) * t224 / 0.2e1 + ((t201 * t268 + t204 * t264) * t248 + (t200 * t268 + t203 * t264) * t247 + (t268 * t233 + t264 * t236 + Icges(2,3)) * t255) * t255 / 0.2e1 + ((-t234 * t265 + t237 * t269 + Icges(1,4)) * V_base(5) + (-t265 * t235 + t269 * t238 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t269 * t234 + t265 * t237 + Icges(1,2)) * V_base(5) + (t235 * t269 + t238 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t255 * (Icges(2,5) * t269 - Icges(2,6) * t265) + V_base(5) * t255 * (Icges(2,5) * t265 + Icges(2,6) * t269) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
