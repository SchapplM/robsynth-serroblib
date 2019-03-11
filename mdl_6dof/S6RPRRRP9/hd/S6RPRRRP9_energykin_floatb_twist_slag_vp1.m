% Calculate kinetic energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:04
% EndTime: 2019-03-09 06:26:07
% DurationCPUTime: 3.61s
% Computational Cost: add. (1469->315), mult. (2139->443), div. (0->0), fcn. (2033->8), ass. (0->154)
t329 = Icges(2,4) + Icges(3,6);
t328 = Icges(2,1) + Icges(3,2);
t327 = Icges(6,1) + Icges(7,1);
t326 = -Icges(3,4) + Icges(2,5);
t325 = Icges(6,4) + Icges(7,4);
t324 = Icges(3,5) - Icges(2,6);
t323 = Icges(7,5) + Icges(6,5);
t322 = Icges(2,2) + Icges(3,3);
t321 = Icges(6,2) + Icges(7,2);
t320 = Icges(7,6) + Icges(6,6);
t319 = Icges(7,3) + Icges(6,3);
t241 = cos(qJ(1));
t318 = t329 * t241;
t238 = sin(qJ(1));
t317 = t329 * t238;
t235 = qJ(4) + qJ(5);
t228 = sin(t235);
t229 = cos(t235);
t237 = sin(qJ(3));
t282 = t237 * t238;
t172 = -t228 * t282 + t229 * t241;
t173 = t228 * t241 + t229 * t282;
t240 = cos(qJ(3));
t279 = t238 * t240;
t316 = t320 * t172 + t323 * t173 - t319 * t279;
t281 = t237 * t241;
t174 = t228 * t281 + t229 * t238;
t175 = t228 * t238 - t229 * t281;
t277 = t240 * t241;
t315 = t320 * t174 + t323 * t175 + t319 * t277;
t314 = t321 * t172 + t325 * t173 - t320 * t279;
t313 = t321 * t174 + t325 * t175 + t320 * t277;
t312 = t325 * t172 + t327 * t173 - t323 * t279;
t311 = t325 * t174 + t327 * t175 + t323 * t277;
t310 = (-t320 * t228 + t323 * t229) * t240 + t319 * t237;
t309 = (-t321 * t228 + t325 * t229) * t240 + t320 * t237;
t308 = (-t325 * t228 + t327 * t229) * t240 + t323 * t237;
t307 = -t322 * t241 - t317;
t306 = t322 * t238 - t318;
t305 = t328 * t238 + t318;
t304 = t328 * t241 - t317;
t273 = pkin(5) * t229;
t301 = -qJ(6) * t240 + t237 * t273;
t239 = cos(qJ(4));
t291 = t239 * pkin(4);
t300 = -pkin(9) * t240 + t237 * t291;
t288 = Icges(4,4) * t237;
t256 = Icges(4,2) * t240 + t288;
t164 = Icges(4,6) * t241 + t238 * t256;
t165 = Icges(4,6) * t238 - t241 * t256;
t287 = Icges(4,4) * t240;
t257 = Icges(4,1) * t237 + t287;
t167 = Icges(4,5) * t241 + t238 * t257;
t168 = Icges(4,5) * t238 - t241 * t257;
t200 = -Icges(4,2) * t237 + t287;
t205 = Icges(4,1) * t240 - t288;
t218 = qJD(3) * t238 + V_base(5);
t219 = qJD(3) * t241 + V_base(4);
t223 = V_base(6) + qJD(1);
t299 = (t164 * t240 + t167 * t237) * t219 + (t165 * t240 + t168 * t237) * t218 + (t200 * t240 + t205 * t237) * t223;
t293 = pkin(7) * t238;
t292 = pkin(7) * t241;
t236 = sin(qJ(4));
t284 = t236 * t238;
t283 = t236 * t241;
t280 = t238 * t239;
t278 = t239 * t241;
t263 = pkin(5) * t228;
t276 = rSges(7,1) * t173 + rSges(7,2) * t172 - rSges(7,3) * t279 + t238 * t301 + t263 * t241;
t275 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t277 + t263 * t238 - t241 * t301;
t274 = (rSges(7,1) * t229 - rSges(7,2) * t228 + t273) * t240 + (qJ(6) + rSges(7,3)) * t237;
t271 = qJD(2) * t241;
t270 = qJD(4) * t240;
t269 = qJD(6) * t240;
t213 = pkin(1) * t241 + qJ(2) * t238;
t268 = t223 * t213 + V_base(2);
t209 = pkin(1) * t238 - qJ(2) * t241;
t267 = V_base(4) * t209 + V_base(3);
t266 = V_base(5) * pkin(6) + V_base(1);
t262 = -t209 - t293;
t261 = qJD(2) * t238 + t266;
t177 = t241 * t270 + t218;
t208 = qJD(4) * t237 + t223;
t260 = V_base(5) * pkin(2) + t261;
t259 = pkin(3) * t237 - pkin(8) * t240;
t258 = rSges(4,1) * t237 + rSges(4,2) * t240;
t255 = Icges(4,5) * t237 + Icges(4,6) * t240;
t251 = (Icges(4,3) * t241 + t238 * t255) * t219 + (Icges(4,3) * t238 - t241 * t255) * t218 + (Icges(4,5) * t240 - Icges(4,6) * t237) * t223;
t250 = t223 * t292 + (-pkin(2) - pkin(6)) * V_base(4) + t268;
t249 = V_base(4) * t293 + (-t213 - t292) * V_base(5) + t267;
t184 = t259 * t238;
t216 = t240 * pkin(3) + t237 * pkin(8);
t248 = t223 * t184 - t216 * t219 + t250;
t185 = t259 * t241;
t247 = t218 * t216 + (t185 + t262) * t223 + t260;
t246 = -t184 * t218 - t219 * t185 + t249;
t140 = pkin(4) * t283 + t238 * t300;
t147 = pkin(9) * t237 + t240 * t291;
t178 = -t238 * t270 + t219;
t245 = t208 * t140 - t147 * t178 + t248;
t141 = pkin(4) * t284 - t241 * t300;
t244 = -t141 * t208 + t177 * t147 + t247;
t243 = -t140 * t177 + t178 * t141 + t246;
t215 = rSges(2,1) * t241 - rSges(2,2) * t238;
t214 = -rSges(3,2) * t241 + rSges(3,3) * t238;
t212 = rSges(4,1) * t240 - rSges(4,2) * t237;
t211 = rSges(2,1) * t238 + rSges(2,2) * t241;
t210 = -rSges(3,2) * t238 - rSges(3,3) * t241;
t191 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t190 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t189 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t186 = qJD(5) * t237 + t208;
t182 = -t237 * t278 + t284;
t181 = t236 * t281 + t280;
t180 = t237 * t280 + t283;
t179 = -t236 * t282 + t278;
t171 = rSges(4,3) * t238 - t241 * t258;
t170 = rSges(5,3) * t237 + (rSges(5,1) * t239 - rSges(5,2) * t236) * t240;
t169 = rSges(4,3) * t241 + t238 * t258;
t166 = Icges(5,5) * t237 + (Icges(5,1) * t239 - Icges(5,4) * t236) * t240;
t163 = Icges(5,6) * t237 + (Icges(5,4) * t239 - Icges(5,2) * t236) * t240;
t160 = Icges(5,3) * t237 + (Icges(5,5) * t239 - Icges(5,6) * t236) * t240;
t158 = (-qJD(4) - qJD(5)) * t279 + t219;
t157 = qJD(5) * t277 + t177;
t156 = rSges(6,3) * t237 + (rSges(6,1) * t229 - rSges(6,2) * t228) * t240;
t146 = V_base(5) * rSges(2,3) - t211 * t223 + t266;
t145 = t215 * t223 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t144 = t211 * V_base(4) - t215 * V_base(5) + V_base(3);
t139 = rSges(5,1) * t182 + rSges(5,2) * t181 + rSges(5,3) * t277;
t138 = rSges(5,1) * t180 + rSges(5,2) * t179 - rSges(5,3) * t279;
t137 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t277;
t136 = Icges(5,1) * t180 + Icges(5,4) * t179 - Icges(5,5) * t279;
t135 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t277;
t134 = Icges(5,4) * t180 + Icges(5,2) * t179 - Icges(5,6) * t279;
t133 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t277;
t132 = Icges(5,5) * t180 + Icges(5,6) * t179 - Icges(5,3) * t279;
t131 = V_base(5) * rSges(3,1) + (-t209 - t210) * t223 + t261;
t130 = -t271 + t214 * t223 + (-rSges(3,1) - pkin(6)) * V_base(4) + t268;
t129 = rSges(6,1) * t175 + rSges(6,2) * t174 + rSges(6,3) * t277;
t127 = rSges(6,1) * t173 + rSges(6,2) * t172 - rSges(6,3) * t279;
t112 = t210 * V_base(4) + (-t213 - t214) * V_base(5) + t267;
t108 = t212 * t218 + (-t171 + t262) * t223 + t260;
t107 = t169 * t223 - t212 * t219 + t250 - t271;
t106 = -t169 * t218 + t171 * t219 + t249;
t105 = -t139 * t208 + t170 * t177 + t247;
t104 = t138 * t208 - t170 * t178 + t248 - t271;
t103 = -t138 * t177 + t139 * t178 + t246;
t102 = -t129 * t186 + t156 * t157 + t244;
t101 = t127 * t186 - t156 * t158 + t245 - t271;
t100 = -t127 * t157 + t129 * t158 + t243;
t99 = t157 * t274 - t186 * t275 - t238 * t269 + t244;
t98 = (-qJD(2) + t269) * t241 + t276 * t186 - t274 * t158 + t245;
t97 = qJD(6) * t237 - t157 * t276 + t158 * t275 + t243;
t1 = t219 * (t299 * t238 + t251 * t241) / 0.2e1 + t218 * (t251 * t238 - t299 * t241) / 0.2e1 + t208 * ((t132 * t178 + t133 * t177 + t160 * t208) * t237 + ((-t134 * t236 + t136 * t239) * t178 + (-t135 * t236 + t137 * t239) * t177 + (-t163 * t236 + t166 * t239) * t208) * t240) / 0.2e1 + t178 * ((-t132 * t279 + t179 * t134 + t180 * t136) * t178 + (-t133 * t279 + t135 * t179 + t137 * t180) * t177 + (-t160 * t279 + t163 * t179 + t166 * t180) * t208) / 0.2e1 + t177 * ((t132 * t277 + t134 * t181 + t136 * t182) * t178 + (t133 * t277 + t181 * t135 + t182 * t137) * t177 + (t160 * t277 + t163 * t181 + t166 * t182) * t208) / 0.2e1 + m(1) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + m(2) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t112 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(4) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(7) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + ((t174 * t309 + t175 * t308 + t277 * t310) * t186 + (t314 * t174 + t312 * t175 + t277 * t316) * t158 + (t313 * t174 + t311 * t175 + t315 * t277) * t157) * t157 / 0.2e1 + ((t172 * t309 + t173 * t308 - t279 * t310) * t186 + (t314 * t172 + t312 * t173 - t316 * t279) * t158 + (t172 * t313 + t173 * t311 - t279 * t315) * t157) * t158 / 0.2e1 + (((-t228 * t309 + t229 * t308) * t186 + (-t228 * t314 + t229 * t312) * t158 + (-t228 * t313 + t229 * t311) * t157) * t240 + (t315 * t157 + t158 * t316 + t310 * t186) * t237) * t186 / 0.2e1 + ((-t164 * t237 + t167 * t240) * t219 + (-t165 * t237 + t168 * t240) * t218 + (-t200 * t237 + t205 * t240 + Icges(3,1) + Icges(2,3)) * t223) * t223 / 0.2e1 + ((t238 * t307 + t241 * t305 + Icges(1,4)) * V_base(5) + (t306 * t238 + t304 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t305 * t238 - t307 * t241 + Icges(1,2)) * V_base(5) + (t238 * t304 - t241 * t306 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t223 * (t324 * t238 + t326 * t241) + V_base(5) * t223 * (t326 * t238 - t324 * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
