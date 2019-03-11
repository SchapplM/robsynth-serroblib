% Calculate kinetic energy for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:08
% EndTime: 2019-03-09 15:21:11
% DurationCPUTime: 3.17s
% Computational Cost: add. (2274->361), mult. (2012->512), div. (0->0), fcn. (1850->12), ass. (0->182)
t339 = Icges(4,3) + Icges(5,3);
t258 = qJ(2) + qJ(3);
t246 = pkin(10) + t258;
t240 = sin(t246);
t241 = cos(t246);
t251 = sin(t258);
t252 = cos(t258);
t338 = Icges(4,5) * t252 + Icges(5,5) * t241 - Icges(4,6) * t251 - Icges(5,6) * t240;
t263 = sin(qJ(1));
t265 = cos(qJ(1));
t321 = Icges(5,4) * t241;
t287 = -Icges(5,2) * t240 + t321;
t164 = -Icges(5,6) * t265 + t263 * t287;
t165 = Icges(5,6) * t263 + t265 * t287;
t322 = Icges(5,4) * t240;
t290 = Icges(5,1) * t241 - t322;
t166 = -Icges(5,5) * t265 + t263 * t290;
t167 = Icges(5,5) * t263 + t265 * t290;
t323 = Icges(4,4) * t252;
t288 = -Icges(4,2) * t251 + t323;
t178 = -Icges(4,6) * t265 + t263 * t288;
t179 = Icges(4,6) * t263 + t265 * t288;
t324 = Icges(4,4) * t251;
t291 = Icges(4,1) * t252 - t324;
t180 = -Icges(4,5) * t265 + t263 * t291;
t181 = Icges(4,5) * t263 + t265 * t291;
t204 = Icges(5,2) * t241 + t322;
t205 = Icges(5,1) * t240 + t321;
t210 = Icges(4,2) * t252 + t324;
t211 = Icges(4,1) * t251 + t323;
t214 = V_base(5) + (-qJD(2) - qJD(3)) * t265;
t239 = qJD(2) * t263 + V_base(4);
t215 = qJD(3) * t263 + t239;
t247 = V_base(6) + qJD(1);
t337 = (-t204 * t240 + t205 * t241 - t210 * t251 + t211 * t252) * t247 + (-t165 * t240 + t167 * t241 - t179 * t251 + t181 * t252) * t215 + (-t164 * t240 + t166 * t241 - t178 * t251 + t180 * t252) * t214;
t336 = (Icges(4,5) * t251 + Icges(5,5) * t240 + Icges(4,6) * t252 + Icges(5,6) * t241) * t247 + (t339 * t263 + t338 * t265) * t215 + (t338 * t263 - t339 * t265) * t214;
t262 = sin(qJ(2));
t332 = pkin(2) * t262;
t331 = pkin(3) * t251;
t264 = cos(qJ(2));
t330 = t264 * pkin(2);
t260 = cos(pkin(11));
t329 = pkin(5) * t260;
t327 = Icges(2,4) * t263;
t326 = Icges(3,4) * t262;
t325 = Icges(3,4) * t264;
t320 = t240 * t263;
t319 = t240 * t265;
t257 = pkin(11) + qJ(6);
t244 = sin(t257);
t318 = t244 * t263;
t317 = t244 * t265;
t245 = cos(t257);
t316 = t245 * t263;
t315 = t245 * t265;
t259 = sin(pkin(11));
t314 = t259 * t263;
t313 = t259 * t265;
t312 = t260 * t263;
t311 = t260 * t265;
t307 = pkin(3) * t252;
t147 = qJ(4) * t263 + t265 * t307;
t293 = pkin(4) * t241 + qJ(5) * t240;
t189 = t293 * t265;
t309 = -t147 - t189;
t172 = -pkin(8) * t265 + t263 * t330;
t236 = t263 * pkin(1) - t265 * pkin(7);
t308 = -t172 - t236;
t305 = qJD(5) * t240;
t304 = qJD(6) * t240;
t303 = V_base(5) * pkin(6) + V_base(1);
t146 = -qJ(4) * t265 + t263 * t307;
t300 = -t146 + t308;
t206 = pkin(4) * t240 - qJ(5) * t241;
t299 = -t206 - t331;
t238 = -qJD(2) * t265 + V_base(5);
t298 = t238 * t332 + t303;
t188 = t293 * t263;
t297 = -t188 + t300;
t296 = rSges(3,1) * t264 - rSges(3,2) * t262;
t295 = rSges(4,1) * t252 - rSges(4,2) * t251;
t294 = rSges(5,1) * t241 - rSges(5,2) * t240;
t292 = Icges(3,1) * t264 - t326;
t289 = -Icges(3,2) * t262 + t325;
t286 = Icges(3,5) * t264 - Icges(3,6) * t262;
t283 = qJD(4) * t263 + t214 * t331 + t298;
t237 = t265 * pkin(1) + t263 * pkin(7);
t282 = -V_base(4) * pkin(6) + t247 * t237 + V_base(2);
t281 = V_base(4) * t236 - t237 * V_base(5) + V_base(3);
t280 = t214 * t206 + t265 * t305 + t283;
t277 = (-Icges(3,3) * t265 + t263 * t286) * t238 + (Icges(3,3) * t263 + t265 * t286) * t239 + (Icges(3,5) * t262 + Icges(3,6) * t264) * t247;
t276 = pkin(9) * t240 + t241 * t329;
t173 = pkin(8) * t263 + t265 * t330;
t275 = t239 * t172 - t173 * t238 + t281;
t274 = t247 * t173 - t239 * t332 + t282;
t273 = t215 * t146 + t275;
t272 = -qJD(4) * t265 + t247 * t147 + t274;
t271 = -qJD(5) * t241 + t215 * t188 + t273;
t270 = t247 * t189 + t263 * t305 + t272;
t192 = -Icges(3,6) * t265 + t263 * t289;
t193 = Icges(3,6) * t263 + t265 * t289;
t194 = -Icges(3,5) * t265 + t263 * t292;
t195 = Icges(3,5) * t263 + t265 * t292;
t227 = Icges(3,2) * t264 + t326;
t230 = Icges(3,1) * t262 + t325;
t267 = (-t193 * t262 + t195 * t264) * t239 + (-t192 * t262 + t194 * t264) * t238 + (-t227 * t262 + t230 * t264) * t247;
t253 = Icges(2,4) * t265;
t235 = rSges(2,1) * t265 - rSges(2,2) * t263;
t234 = rSges(2,1) * t263 + rSges(2,2) * t265;
t233 = rSges(3,1) * t262 + rSges(3,2) * t264;
t232 = Icges(2,1) * t265 - t327;
t231 = Icges(2,1) * t263 + t253;
t229 = -Icges(2,2) * t263 + t253;
t228 = Icges(2,2) * t265 + t327;
t221 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t220 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t219 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = -qJD(6) * t241 + t247;
t212 = rSges(4,1) * t251 + rSges(4,2) * t252;
t207 = rSges(5,1) * t240 + rSges(5,2) * t241;
t201 = rSges(3,3) * t263 + t265 * t296;
t200 = -rSges(3,3) * t265 + t263 * t296;
t199 = t241 * t311 + t314;
t198 = -t241 * t313 + t312;
t197 = t241 * t312 - t313;
t196 = -t241 * t314 - t311;
t187 = rSges(4,3) * t263 + t265 * t295;
t186 = -rSges(4,3) * t265 + t263 * t295;
t185 = t241 * t315 + t318;
t184 = -t241 * t317 + t316;
t183 = t241 * t316 - t317;
t182 = -t241 * t318 - t315;
t175 = t265 * t304 + t215;
t174 = t263 * t304 + t214;
t171 = rSges(5,3) * t263 + t265 * t294;
t170 = -rSges(5,3) * t265 + t263 * t294;
t169 = V_base(5) * rSges(2,3) - t234 * t247 + t303;
t168 = t235 * t247 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t160 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t157 = -rSges(6,3) * t241 + (rSges(6,1) * t260 - rSges(6,2) * t259) * t240;
t156 = -Icges(6,5) * t241 + (Icges(6,1) * t260 - Icges(6,4) * t259) * t240;
t155 = -Icges(6,6) * t241 + (Icges(6,4) * t260 - Icges(6,2) * t259) * t240;
t154 = -Icges(6,3) * t241 + (Icges(6,5) * t260 - Icges(6,6) * t259) * t240;
t152 = -rSges(7,3) * t241 + (rSges(7,1) * t245 - rSges(7,2) * t244) * t240;
t151 = -Icges(7,5) * t241 + (Icges(7,1) * t245 - Icges(7,4) * t244) * t240;
t150 = -Icges(7,6) * t241 + (Icges(7,4) * t245 - Icges(7,2) * t244) * t240;
t149 = -Icges(7,3) * t241 + (Icges(7,5) * t245 - Icges(7,6) * t244) * t240;
t145 = -pkin(9) * t241 + t240 * t329;
t142 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t319;
t141 = rSges(6,1) * t197 + rSges(6,2) * t196 + rSges(6,3) * t320;
t140 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t319;
t139 = Icges(6,1) * t197 + Icges(6,4) * t196 + Icges(6,5) * t320;
t138 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t319;
t137 = Icges(6,4) * t197 + Icges(6,2) * t196 + Icges(6,6) * t320;
t136 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t319;
t135 = Icges(6,5) * t197 + Icges(6,6) * t196 + Icges(6,3) * t320;
t134 = pkin(5) * t314 + t265 * t276;
t133 = -pkin(5) * t313 + t263 * t276;
t132 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t319;
t131 = rSges(7,1) * t183 + rSges(7,2) * t182 + rSges(7,3) * t320;
t130 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t319;
t129 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t320;
t128 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t319;
t127 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t320;
t126 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t319;
t125 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t320;
t124 = t233 * t238 + (-t200 - t236) * t247 + t303;
t123 = t201 * t247 - t233 * t239 + t282;
t122 = t200 * t239 - t201 * t238 + t281;
t121 = t212 * t214 + (-t186 + t308) * t247 + t298;
t120 = t187 * t247 - t212 * t215 + t274;
t119 = t186 * t215 - t187 * t214 + t275;
t118 = t207 * t214 + (-t170 + t300) * t247 + t283;
t117 = t171 * t247 + (-t207 - t331) * t215 + t272;
t116 = t170 * t215 + (-t147 - t171) * t214 + t273;
t115 = t157 * t214 + (-t141 + t297) * t247 + t280;
t114 = t142 * t247 + (-t157 + t299) * t215 + t270;
t113 = t141 * t215 + (-t142 + t309) * t214 + t271;
t112 = -t131 * t213 + t145 * t214 + t152 * t174 + (-t133 + t297) * t247 + t280;
t111 = t132 * t213 + t134 * t247 - t152 * t175 + (-t145 + t299) * t215 + t270;
t110 = t131 * t175 - t132 * t174 + t133 * t215 + (-t134 + t309) * t214 + t271;
t1 = t175 * ((t126 * t319 + t184 * t128 + t185 * t130) * t175 + (t125 * t319 + t127 * t184 + t129 * t185) * t174 + (t149 * t319 + t150 * t184 + t151 * t185) * t213) / 0.2e1 + t174 * ((t126 * t320 + t128 * t182 + t130 * t183) * t175 + (t125 * t320 + t182 * t127 + t183 * t129) * t174 + (t149 * t320 + t150 * t182 + t151 * t183) * t213) / 0.2e1 + t239 * (t263 * t277 + t265 * t267) / 0.2e1 + t238 * (t263 * t267 - t265 * t277) / 0.2e1 + m(2) * (t160 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t213 * ((-t125 * t174 - t126 * t175 - t149 * t213) * t241 + ((-t128 * t244 + t130 * t245) * t175 + (-t127 * t244 + t129 * t245) * t174 + (-t150 * t244 + t151 * t245) * t213) * t240) / 0.2e1 + m(1) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + ((-t228 * t263 + t231 * t265 + Icges(1,4)) * V_base(5) + (-t263 * t229 + t265 * t232 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t265 * t228 + t263 * t231 + Icges(1,2)) * V_base(5) + (t229 * t265 + t232 * t263 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t320 + t138 * t196 + t140 * t197) * t215 + (t135 * t320 + t196 * t137 + t197 * t139) * t214 + (t154 * t320 + t155 * t196 + t156 * t197) * t247 - t336 * t265 + t337 * t263) * t214 / 0.2e1 + ((t136 * t319 + t198 * t138 + t199 * t140) * t215 + (t135 * t319 + t137 * t198 + t139 * t199) * t214 + (t154 * t319 + t155 * t198 + t156 * t199) * t247 + t337 * t265 + t336 * t263) * t215 / 0.2e1 + ((t193 * t264 + t195 * t262) * t239 + (t192 * t264 + t194 * t262) * t238 + (t179 * t252 + t181 * t251 + (-t136 + t165) * t241 + (-t138 * t259 + t140 * t260 + t167) * t240) * t215 + (t178 * t252 + t180 * t251 + (-t135 + t164) * t241 + (-t137 * t259 + t139 * t260 + t166) * t240) * t214 + (t252 * t210 + t251 * t211 + t264 * t227 + t262 * t230 + Icges(2,3) + (-t154 + t204) * t241 + (-t155 * t259 + t156 * t260 + t205) * t240) * t247) * t247 / 0.2e1 + t247 * V_base(4) * (Icges(2,5) * t265 - Icges(2,6) * t263) + t247 * V_base(5) * (Icges(2,5) * t263 + Icges(2,6) * t265) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
