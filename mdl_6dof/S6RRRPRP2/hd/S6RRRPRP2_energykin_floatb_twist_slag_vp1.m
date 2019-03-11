% Calculate kinetic energy for
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:02
% EndTime: 2019-03-09 16:35:05
% DurationCPUTime: 2.75s
% Computational Cost: add. (2155->306), mult. (2038->432), div. (0->0), fcn. (1886->10), ass. (0->161)
t337 = Icges(6,1) + Icges(7,1);
t336 = -Icges(6,4) + Icges(7,5);
t335 = Icges(7,4) + Icges(6,5);
t334 = Icges(6,2) + Icges(7,3);
t333 = -Icges(7,6) + Icges(6,6);
t332 = Icges(4,3) + Icges(5,3);
t331 = -Icges(6,3) - Icges(7,2);
t245 = qJ(2) + qJ(3);
t234 = pkin(10) + t245;
t230 = sin(t234);
t231 = cos(t234);
t239 = sin(t245);
t240 = cos(t245);
t330 = Icges(4,5) * t240 + Icges(5,5) * t231 - Icges(4,6) * t239 - Icges(5,6) * t230;
t329 = rSges(7,1) + pkin(5);
t328 = rSges(7,3) + qJ(6);
t249 = cos(qJ(5));
t251 = cos(qJ(1));
t294 = t249 * t251;
t246 = sin(qJ(5));
t248 = sin(qJ(1));
t297 = t246 * t248;
t191 = t231 * t297 + t294;
t295 = t248 * t249;
t296 = t246 * t251;
t192 = t231 * t295 - t296;
t299 = t230 * t248;
t327 = t334 * t191 + t336 * t192 - t333 * t299;
t193 = t231 * t296 - t295;
t194 = t231 * t294 + t297;
t298 = t230 * t251;
t326 = t334 * t193 + t336 * t194 - t333 * t298;
t325 = -t333 * t191 + t335 * t192 - t331 * t299;
t324 = -t333 * t193 + t335 * t194 - t331 * t298;
t323 = t336 * t191 + t337 * t192 + t335 * t299;
t322 = t336 * t193 + t337 * t194 + t335 * t298;
t321 = t333 * t231 + (t334 * t246 + t336 * t249) * t230;
t320 = t331 * t231 + (-t333 * t246 + t335 * t249) * t230;
t319 = -t335 * t231 + (t336 * t246 + t337 * t249) * t230;
t300 = Icges(5,4) * t231;
t272 = -Icges(5,2) * t230 + t300;
t159 = -Icges(5,6) * t251 + t248 * t272;
t160 = Icges(5,6) * t248 + t251 * t272;
t301 = Icges(5,4) * t230;
t275 = Icges(5,1) * t231 - t301;
t161 = -Icges(5,5) * t251 + t248 * t275;
t162 = Icges(5,5) * t248 + t251 * t275;
t302 = Icges(4,4) * t240;
t273 = -Icges(4,2) * t239 + t302;
t173 = -Icges(4,6) * t251 + t248 * t273;
t174 = Icges(4,6) * t248 + t251 * t273;
t303 = Icges(4,4) * t239;
t276 = Icges(4,1) * t240 - t303;
t175 = -Icges(4,5) * t251 + t248 * t276;
t176 = Icges(4,5) * t248 + t251 * t276;
t196 = Icges(5,2) * t231 + t301;
t197 = Icges(5,1) * t230 + t300;
t202 = Icges(4,2) * t240 + t303;
t203 = Icges(4,1) * t239 + t302;
t206 = V_base(5) + (-qJD(2) - qJD(3)) * t251;
t229 = qJD(2) * t248 + V_base(4);
t207 = qJD(3) * t248 + t229;
t235 = V_base(6) + qJD(1);
t318 = (-t196 * t230 + t197 * t231 - t202 * t239 + t203 * t240) * t235 + (-t160 * t230 + t162 * t231 - t174 * t239 + t176 * t240) * t207 + (-t159 * t230 + t161 * t231 - t173 * t239 + t175 * t240) * t206;
t317 = (Icges(4,5) * t239 + Icges(5,5) * t230 + Icges(4,6) * t240 + Icges(5,6) * t231) * t235 + (t332 * t248 + t330 * t251) * t207 + (t330 * t248 - t332 * t251) * t206;
t247 = sin(qJ(2));
t310 = pkin(2) * t247;
t309 = pkin(3) * t239;
t250 = cos(qJ(2));
t308 = t250 * pkin(2);
t306 = Icges(2,4) * t248;
t305 = Icges(3,4) * t247;
t304 = Icges(3,4) * t250;
t293 = rSges(7,2) * t299 + t328 * t191 + t192 * t329;
t292 = rSges(7,2) * t298 + t328 * t193 + t194 * t329;
t291 = -rSges(7,2) * t231 + (t328 * t246 + t249 * t329) * t230;
t167 = -pkin(8) * t251 + t248 * t308;
t226 = t248 * pkin(1) - t251 * pkin(7);
t290 = -t167 - t226;
t289 = pkin(3) * t240;
t287 = qJD(5) * t230;
t286 = V_base(5) * pkin(6) + V_base(1);
t141 = -qJ(4) * t251 + t248 * t289;
t283 = -t141 + t290;
t228 = -qJD(2) * t251 + V_base(5);
t282 = t228 * t310 + t286;
t281 = pkin(4) * t231 + pkin(9) * t230;
t280 = rSges(3,1) * t250 - rSges(3,2) * t247;
t279 = rSges(4,1) * t240 - rSges(4,2) * t239;
t278 = rSges(5,1) * t231 - rSges(5,2) * t230;
t277 = Icges(3,1) * t250 - t305;
t274 = -Icges(3,2) * t247 + t304;
t271 = Icges(3,5) * t250 - Icges(3,6) * t247;
t268 = qJD(4) * t248 + t206 * t309 + t282;
t227 = t251 * pkin(1) + t248 * pkin(7);
t267 = -V_base(4) * pkin(6) + t235 * t227 + V_base(2);
t266 = V_base(4) * t226 - t227 * V_base(5) + V_base(3);
t263 = (-Icges(3,3) * t251 + t248 * t271) * t228 + (Icges(3,3) * t248 + t251 * t271) * t229 + (Icges(3,5) * t247 + Icges(3,6) * t250) * t235;
t168 = pkin(8) * t248 + t251 * t308;
t262 = t229 * t167 - t168 * t228 + t266;
t261 = t235 * t168 - t229 * t310 + t267;
t260 = t207 * t141 + t262;
t180 = t281 * t248;
t199 = pkin(4) * t230 - pkin(9) * t231;
t259 = t206 * t199 + (-t180 + t283) * t235 + t268;
t142 = qJ(4) * t248 + t251 * t289;
t258 = -qJD(4) * t251 + t235 * t142 + t261;
t181 = t281 * t251;
t257 = t207 * t180 + (-t142 - t181) * t206 + t260;
t256 = t235 * t181 + (-t199 - t309) * t207 + t258;
t184 = -Icges(3,6) * t251 + t248 * t274;
t185 = Icges(3,6) * t248 + t251 * t274;
t186 = -Icges(3,5) * t251 + t248 * t277;
t187 = Icges(3,5) * t248 + t251 * t277;
t217 = Icges(3,2) * t250 + t305;
t220 = Icges(3,1) * t247 + t304;
t253 = (-t185 * t247 + t187 * t250) * t229 + (-t184 * t247 + t186 * t250) * t228 + (-t217 * t247 + t220 * t250) * t235;
t241 = Icges(2,4) * t251;
t225 = rSges(2,1) * t251 - rSges(2,2) * t248;
t224 = rSges(2,1) * t248 + rSges(2,2) * t251;
t223 = rSges(3,1) * t247 + rSges(3,2) * t250;
t222 = Icges(2,1) * t251 - t306;
t221 = Icges(2,1) * t248 + t241;
t219 = -Icges(2,2) * t248 + t241;
t218 = Icges(2,2) * t251 + t306;
t213 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t212 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t211 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t205 = -qJD(5) * t231 + t235;
t204 = rSges(4,1) * t239 + rSges(4,2) * t240;
t198 = rSges(5,1) * t230 + rSges(5,2) * t231;
t189 = rSges(3,3) * t248 + t251 * t280;
t188 = -rSges(3,3) * t251 + t248 * t280;
t178 = rSges(4,3) * t248 + t251 * t279;
t177 = -rSges(4,3) * t251 + t248 * t279;
t170 = t251 * t287 + t207;
t169 = t248 * t287 + t206;
t166 = rSges(5,3) * t248 + t251 * t278;
t165 = -rSges(5,3) * t251 + t248 * t278;
t164 = V_base(5) * rSges(2,3) - t224 * t235 + t286;
t163 = t225 * t235 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t224 * V_base(4) - t225 * V_base(5) + V_base(3);
t152 = -rSges(6,3) * t231 + (rSges(6,1) * t249 - rSges(6,2) * t246) * t230;
t136 = rSges(6,1) * t194 - rSges(6,2) * t193 + rSges(6,3) * t298;
t134 = rSges(6,1) * t192 - rSges(6,2) * t191 + rSges(6,3) * t299;
t120 = t223 * t228 + (-t188 - t226) * t235 + t286;
t119 = t189 * t235 - t223 * t229 + t267;
t118 = t188 * t229 - t189 * t228 + t266;
t117 = t204 * t206 + (-t177 + t290) * t235 + t282;
t116 = t178 * t235 - t204 * t207 + t261;
t115 = t177 * t207 - t178 * t206 + t262;
t114 = t198 * t206 + (-t165 + t283) * t235 + t268;
t113 = t166 * t235 + (-t198 - t309) * t207 + t258;
t112 = t165 * t207 + (-t142 - t166) * t206 + t260;
t111 = -t134 * t205 + t152 * t169 + t259;
t110 = t136 * t205 - t152 * t170 + t256;
t109 = t134 * t170 - t136 * t169 + t257;
t108 = qJD(6) * t193 + t169 * t291 - t205 * t293 + t259;
t107 = qJD(6) * t191 - t170 * t291 + t205 * t292 + t256;
t106 = qJD(6) * t230 * t246 - t169 * t292 + t170 * t293 + t257;
t1 = m(1) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(3) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t229 * (t248 * t263 + t251 * t253) / 0.2e1 + t228 * (t248 * t253 - t251 * t263) / 0.2e1 + ((t191 * t321 + t192 * t319 + t299 * t320) * t205 + (t191 * t326 + t192 * t322 + t299 * t324) * t170 + (t327 * t191 + t323 * t192 + t325 * t299) * t169) * t169 / 0.2e1 + ((t193 * t321 + t194 * t319 + t298 * t320) * t205 + (t326 * t193 + t322 * t194 + t324 * t298) * t170 + (t193 * t327 + t323 * t194 + t325 * t298) * t169) * t170 / 0.2e1 + ((-t169 * t325 - t170 * t324 - t320 * t205) * t231 + ((t246 * t321 + t319 * t249) * t205 + (t246 * t326 + t249 * t322) * t170 + (t246 * t327 + t323 * t249) * t169) * t230) * t205 / 0.2e1 + (t318 * t248 - t317 * t251) * t206 / 0.2e1 + (t317 * t248 + t318 * t251) * t207 / 0.2e1 + ((-t218 * t248 + t221 * t251 + Icges(1,4)) * V_base(5) + (-t219 * t248 + t222 * t251 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t218 * t251 + t221 * t248 + Icges(1,2)) * V_base(5) + (t219 * t251 + t222 * t248 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t185 * t250 + t187 * t247) * t229 + (t184 * t250 + t186 * t247) * t228 + (t160 * t231 + t162 * t230 + t174 * t240 + t176 * t239) * t207 + (t159 * t231 + t161 * t230 + t173 * t240 + t175 * t239) * t206 + (t196 * t231 + t197 * t230 + t202 * t240 + t203 * t239 + t217 * t250 + t220 * t247 + Icges(2,3)) * t235) * t235 / 0.2e1 + V_base(4) * t235 * (Icges(2,5) * t251 - Icges(2,6) * t248) + V_base(5) * t235 * (Icges(2,5) * t248 + Icges(2,6) * t251) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
