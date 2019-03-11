% Calculate kinetic energy for
% S6RRRPRP1
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
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:30:39
% EndTime: 2019-03-09 16:30:41
% DurationCPUTime: 2.80s
% Computational Cost: add. (2197->311), mult. (2047->439), div. (0->0), fcn. (1885->10), ass. (0->162)
t338 = Icges(6,1) + Icges(7,1);
t337 = Icges(6,4) + Icges(7,4);
t336 = -Icges(7,5) - Icges(6,5);
t335 = Icges(6,2) + Icges(7,2);
t334 = -Icges(7,6) - Icges(6,6);
t333 = Icges(4,3) + Icges(5,3);
t332 = -Icges(7,3) - Icges(6,3);
t243 = qJ(2) + qJ(3);
t232 = pkin(10) + t243;
t228 = sin(t232);
t229 = cos(t232);
t237 = sin(t243);
t238 = cos(t243);
t331 = Icges(4,5) * t238 + Icges(5,5) * t229 - Icges(4,6) * t237 - Icges(5,6) * t228;
t248 = cos(qJ(5));
t250 = cos(qJ(1));
t295 = t248 * t250;
t245 = sin(qJ(5));
t247 = sin(qJ(1));
t298 = t245 * t247;
t189 = -t229 * t298 - t295;
t296 = t247 * t248;
t297 = t245 * t250;
t190 = t229 * t296 - t297;
t300 = t228 * t247;
t330 = -t189 * t334 - t190 * t336 - t300 * t332;
t191 = -t229 * t297 + t296;
t192 = t229 * t295 + t298;
t299 = t228 * t250;
t329 = -t191 * t334 - t192 * t336 - t299 * t332;
t328 = t189 * t335 + t190 * t337 - t300 * t334;
t327 = t191 * t335 + t192 * t337 - t299 * t334;
t326 = t337 * t189 + t190 * t338 - t336 * t300;
t325 = t337 * t191 + t192 * t338 - t336 * t299;
t324 = t332 * t229 + (t245 * t334 - t248 * t336) * t228;
t323 = t334 * t229 + (-t245 * t335 + t248 * t337) * t228;
t322 = t336 * t229 + (-t337 * t245 + t248 * t338) * t228;
t301 = Icges(5,4) * t229;
t272 = -Icges(5,2) * t228 + t301;
t158 = -Icges(5,6) * t250 + t247 * t272;
t159 = Icges(5,6) * t247 + t250 * t272;
t302 = Icges(5,4) * t228;
t275 = Icges(5,1) * t229 - t302;
t160 = -Icges(5,5) * t250 + t247 * t275;
t161 = Icges(5,5) * t247 + t250 * t275;
t303 = Icges(4,4) * t238;
t273 = -Icges(4,2) * t237 + t303;
t172 = -Icges(4,6) * t250 + t247 * t273;
t173 = Icges(4,6) * t247 + t250 * t273;
t304 = Icges(4,4) * t237;
t276 = Icges(4,1) * t238 - t304;
t174 = -Icges(4,5) * t250 + t247 * t276;
t175 = Icges(4,5) * t247 + t250 * t276;
t194 = Icges(5,2) * t229 + t302;
t195 = Icges(5,1) * t228 + t301;
t200 = Icges(4,2) * t238 + t304;
t201 = Icges(4,1) * t237 + t303;
t204 = V_base(5) + (-qJD(2) - qJD(3)) * t250;
t227 = qJD(2) * t247 + V_base(4);
t205 = qJD(3) * t247 + t227;
t233 = V_base(6) + qJD(1);
t321 = (-t194 * t228 + t195 * t229 - t200 * t237 + t201 * t238) * t233 + (-t159 * t228 + t161 * t229 - t173 * t237 + t175 * t238) * t205 + (-t158 * t228 + t160 * t229 - t172 * t237 + t174 * t238) * t204;
t320 = (Icges(4,5) * t237 + Icges(5,5) * t228 + Icges(4,6) * t238 + Icges(5,6) * t229) * t233 + (t247 * t333 + t250 * t331) * t205 + (t247 * t331 - t250 * t333) * t204;
t246 = sin(qJ(2));
t313 = pkin(2) * t246;
t312 = pkin(3) * t237;
t249 = cos(qJ(2));
t311 = t249 * pkin(2);
t310 = pkin(5) * t248;
t307 = Icges(2,4) * t247;
t306 = Icges(3,4) * t246;
t305 = Icges(3,4) * t249;
t265 = qJ(6) * t228 + t229 * t310;
t294 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t300 - pkin(5) * t297 + t247 * t265;
t293 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t299 + pkin(5) * t298 + t250 * t265;
t292 = (-qJ(6) - rSges(7,3)) * t229 + (rSges(7,1) * t248 - rSges(7,2) * t245 + t310) * t228;
t166 = -pkin(8) * t250 + t247 * t311;
t224 = t247 * pkin(1) - t250 * pkin(7);
t291 = -t166 - t224;
t290 = pkin(3) * t238;
t288 = qJD(5) * t228;
t287 = qJD(6) * t228;
t286 = V_base(5) * pkin(6) + V_base(1);
t140 = -qJ(4) * t250 + t247 * t290;
t283 = -t140 + t291;
t226 = -qJD(2) * t250 + V_base(5);
t282 = t226 * t313 + t286;
t281 = pkin(4) * t229 + pkin(9) * t228;
t280 = rSges(3,1) * t249 - rSges(3,2) * t246;
t279 = rSges(4,1) * t238 - rSges(4,2) * t237;
t278 = rSges(5,1) * t229 - rSges(5,2) * t228;
t277 = Icges(3,1) * t249 - t306;
t274 = -Icges(3,2) * t246 + t305;
t271 = Icges(3,5) * t249 - Icges(3,6) * t246;
t268 = qJD(4) * t247 + t204 * t312 + t282;
t225 = t250 * pkin(1) + t247 * pkin(7);
t267 = -V_base(4) * pkin(6) + t233 * t225 + V_base(2);
t266 = V_base(4) * t224 - t225 * V_base(5) + V_base(3);
t262 = (-Icges(3,3) * t250 + t247 * t271) * t226 + (Icges(3,3) * t247 + t250 * t271) * t227 + (Icges(3,5) * t246 + Icges(3,6) * t249) * t233;
t167 = pkin(8) * t247 + t250 * t311;
t261 = t227 * t166 - t167 * t226 + t266;
t260 = t233 * t167 - t227 * t313 + t267;
t259 = t205 * t140 + t261;
t178 = t281 * t247;
t197 = pkin(4) * t228 - pkin(9) * t229;
t258 = t204 * t197 + (-t178 + t283) * t233 + t268;
t141 = qJ(4) * t247 + t250 * t290;
t257 = -qJD(4) * t250 + t233 * t141 + t260;
t179 = t281 * t250;
t256 = t205 * t178 + (-t141 - t179) * t204 + t259;
t255 = t233 * t179 + (-t197 - t312) * t205 + t257;
t182 = -Icges(3,6) * t250 + t247 * t274;
t183 = Icges(3,6) * t247 + t250 * t274;
t184 = -Icges(3,5) * t250 + t247 * t277;
t185 = Icges(3,5) * t247 + t250 * t277;
t215 = Icges(3,2) * t249 + t306;
t218 = Icges(3,1) * t246 + t305;
t252 = (-t183 * t246 + t185 * t249) * t227 + (-t182 * t246 + t184 * t249) * t226 + (-t215 * t246 + t218 * t249) * t233;
t239 = Icges(2,4) * t250;
t223 = rSges(2,1) * t250 - rSges(2,2) * t247;
t222 = rSges(2,1) * t247 + rSges(2,2) * t250;
t221 = rSges(3,1) * t246 + rSges(3,2) * t249;
t220 = Icges(2,1) * t250 - t307;
t219 = Icges(2,1) * t247 + t239;
t217 = -Icges(2,2) * t247 + t239;
t216 = Icges(2,2) * t250 + t307;
t211 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t210 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t209 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t203 = -qJD(5) * t229 + t233;
t202 = rSges(4,1) * t237 + rSges(4,2) * t238;
t196 = rSges(5,1) * t228 + rSges(5,2) * t229;
t187 = rSges(3,3) * t247 + t250 * t280;
t186 = -rSges(3,3) * t250 + t247 * t280;
t177 = rSges(4,3) * t247 + t250 * t279;
t176 = -rSges(4,3) * t250 + t247 * t279;
t169 = t250 * t288 + t205;
t168 = t247 * t288 + t204;
t165 = rSges(5,3) * t247 + t250 * t278;
t164 = -rSges(5,3) * t250 + t247 * t278;
t163 = V_base(5) * rSges(2,3) - t222 * t233 + t286;
t162 = t223 * t233 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = t222 * V_base(4) - t223 * V_base(5) + V_base(3);
t151 = -rSges(6,3) * t229 + (rSges(6,1) * t248 - rSges(6,2) * t245) * t228;
t136 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t299;
t134 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t300;
t118 = t221 * t226 + (-t186 - t224) * t233 + t286;
t117 = t187 * t233 - t221 * t227 + t267;
t116 = t186 * t227 - t187 * t226 + t266;
t115 = t202 * t204 + (-t176 + t291) * t233 + t282;
t114 = t177 * t233 - t202 * t205 + t260;
t113 = t176 * t205 - t177 * t204 + t261;
t112 = t196 * t204 + (-t164 + t283) * t233 + t268;
t111 = t165 * t233 + (-t196 - t312) * t205 + t257;
t110 = t164 * t205 + (-t141 - t165) * t204 + t259;
t109 = -t134 * t203 + t151 * t168 + t258;
t108 = t136 * t203 - t151 * t169 + t255;
t107 = t134 * t169 - t136 * t168 + t256;
t106 = t168 * t292 - t203 * t294 + t250 * t287 + t258;
t105 = -t169 * t292 + t203 * t293 + t247 * t287 + t255;
t104 = -qJD(6) * t229 - t168 * t293 + t169 * t294 + t256;
t1 = m(1) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(2) * (t154 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t227 * (t247 * t262 + t250 * t252) / 0.2e1 + t226 * (t247 * t252 - t250 * t262) / 0.2e1 + ((t189 * t323 + t190 * t322 + t300 * t324) * t203 + (t189 * t327 + t190 * t325 + t300 * t329) * t169 + (t328 * t189 + t326 * t190 + t330 * t300) * t168) * t168 / 0.2e1 + ((t191 * t323 + t192 * t322 + t299 * t324) * t203 + (t327 * t191 + t325 * t192 + t329 * t299) * t169 + (t328 * t191 + t326 * t192 + t330 * t299) * t168) * t169 / 0.2e1 + ((-t168 * t330 - t329 * t169 - t324 * t203) * t229 + ((-t245 * t323 + t248 * t322) * t203 + (-t245 * t327 + t248 * t325) * t169 + (-t245 * t328 + t248 * t326) * t168) * t228) * t203 / 0.2e1 + (t321 * t247 - t320 * t250) * t204 / 0.2e1 + (t320 * t247 + t321 * t250) * t205 / 0.2e1 + ((-t216 * t247 + t219 * t250 + Icges(1,4)) * V_base(5) + (-t217 * t247 + t220 * t250 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t216 * t250 + t219 * t247 + Icges(1,2)) * V_base(5) + (t217 * t250 + t220 * t247 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t183 * t249 + t185 * t246) * t227 + (t182 * t249 + t184 * t246) * t226 + (t159 * t229 + t161 * t228 + t173 * t238 + t175 * t237) * t205 + (t158 * t229 + t160 * t228 + t172 * t238 + t174 * t237) * t204 + (t194 * t229 + t195 * t228 + t200 * t238 + t201 * t237 + t215 * t249 + t218 * t246 + Icges(2,3)) * t233) * t233 / 0.2e1 + t233 * V_base(4) * (Icges(2,5) * t250 - Icges(2,6) * t247) + t233 * V_base(5) * (Icges(2,5) * t247 + Icges(2,6) * t250) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
