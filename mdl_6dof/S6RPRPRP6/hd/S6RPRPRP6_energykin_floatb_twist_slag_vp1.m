% Calculate kinetic energy for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:17:53
% EndTime: 2019-03-09 03:17:57
% DurationCPUTime: 3.82s
% Computational Cost: add. (1598->295), mult. (1918->403), div. (0->0), fcn. (1756->8), ass. (0->156)
t346 = Icges(4,4) + Icges(5,6);
t345 = Icges(4,1) + Icges(5,2);
t344 = -Icges(4,2) - Icges(5,3);
t235 = pkin(9) + qJ(3);
t229 = cos(t235);
t343 = t346 * t229;
t228 = sin(t235);
t342 = t346 * t228;
t341 = Icges(5,4) - Icges(4,5);
t340 = Icges(5,5) - Icges(4,6);
t339 = t344 * t228 + t343;
t338 = t345 * t229 - t342;
t337 = Icges(5,1) + Icges(4,3);
t336 = Icges(6,1) + Icges(7,1);
t335 = Icges(6,4) + Icges(7,4);
t334 = Icges(7,5) + Icges(6,5);
t333 = Icges(6,2) + Icges(7,2);
t332 = Icges(7,6) + Icges(6,6);
t331 = Icges(7,3) + Icges(6,3);
t241 = sin(qJ(1));
t243 = cos(qJ(1));
t330 = t339 * t241 + t340 * t243;
t329 = -t340 * t241 + t339 * t243;
t328 = t338 * t241 + t341 * t243;
t327 = -t341 * t241 + t338 * t243;
t326 = t344 * t229 - t342;
t325 = t345 * t228 + t343;
t324 = t340 * t228 - t341 * t229;
t242 = cos(qJ(5));
t288 = t242 * t243;
t240 = sin(qJ(5));
t290 = t241 * t240;
t185 = t228 * t288 - t290;
t289 = t241 * t242;
t291 = t240 * t243;
t186 = t228 * t291 + t289;
t292 = t229 * t243;
t323 = t332 * t185 + t334 * t186 + t331 * t292;
t187 = t228 * t289 + t291;
t188 = t228 * t290 - t288;
t293 = t229 * t241;
t322 = t332 * t187 + t334 * t188 + t331 * t293;
t321 = t333 * t185 + t335 * t186 + t332 * t292;
t320 = t333 * t187 + t335 * t188 + t332 * t293;
t319 = t335 * t185 + t336 * t186 + t334 * t292;
t318 = t335 * t187 + t336 * t188 + t334 * t293;
t317 = (-t334 * t240 - t332 * t242) * t229 + t331 * t228;
t316 = (-t335 * t240 - t333 * t242) * t229 + t332 * t228;
t315 = (-t336 * t240 - t335 * t242) * t229 + t334 * t228;
t223 = -qJD(3) * t243 + V_base(5);
t224 = qJD(3) * t241 + V_base(4);
t230 = V_base(6) + qJD(1);
t314 = (t326 * t228 + t325 * t229) * t230 + (-t329 * t228 + t327 * t229) * t224 + (-t330 * t228 + t328 * t229) * t223;
t313 = (-t341 * t228 - t340 * t229) * t230 + (t337 * t241 + t324 * t243) * t224 + (t324 * t241 - t337 * t243) * t223;
t236 = sin(pkin(9));
t306 = pkin(2) * t236;
t305 = pkin(5) * t240;
t304 = pkin(8) * t228;
t237 = cos(pkin(9));
t303 = pkin(2) * t237;
t302 = pkin(5) * t242;
t300 = Icges(2,4) * t241;
t299 = Icges(3,4) * t236;
t298 = Icges(3,4) * t237;
t257 = qJ(6) * t229 + t228 * t305;
t286 = t186 * rSges(7,1) + t185 * rSges(7,2) + rSges(7,3) * t292 + t241 * t302 + t243 * t257;
t285 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t293 + t241 * t257 - t243 * t302;
t284 = (-rSges(7,1) * t240 - rSges(7,2) * t242 - t305) * t229 + (rSges(7,3) + qJ(6)) * t228;
t150 = -pkin(7) * t243 + t241 * t303;
t219 = t241 * pkin(1) - qJ(2) * t243;
t283 = -t150 - t219;
t282 = qJD(4) * t228;
t281 = qJD(5) * t229;
t280 = qJD(6) * t229;
t279 = V_base(4) * t219 + V_base(3);
t278 = V_base(5) * pkin(6) + V_base(1);
t269 = pkin(3) * t229 + qJ(4) * t228;
t181 = t269 * t241;
t275 = -t181 + t283;
t274 = qJD(2) * t241 + t278;
t273 = V_base(5) * t306 + t274;
t272 = rSges(3,1) * t237 - rSges(3,2) * t236;
t271 = rSges(4,1) * t229 - rSges(4,2) * t228;
t270 = -rSges(5,2) * t229 + rSges(5,3) * t228;
t268 = Icges(3,1) * t237 - t299;
t266 = -Icges(3,2) * t236 + t298;
t263 = Icges(3,5) * t237 - Icges(3,6) * t236;
t221 = pkin(1) * t243 + t241 * qJ(2);
t259 = -qJD(2) * t243 + t230 * t221 + V_base(2);
t197 = pkin(3) * t228 - qJ(4) * t229;
t258 = t223 * t197 + t243 * t282 + t273;
t151 = pkin(7) * t241 + t243 * t303;
t254 = V_base(4) * t150 + (-t151 - t221) * V_base(5) + t279;
t253 = (-Icges(3,3) * t243 + t241 * t263) * V_base(5) + (Icges(3,3) * t241 + t243 * t263) * V_base(4) + (Icges(3,5) * t236 + Icges(3,6) * t237) * t230;
t252 = t230 * t151 + (-pkin(6) - t306) * V_base(4) + t259;
t251 = -qJD(4) * t229 + t224 * t181 + t254;
t201 = -pkin(4) * t243 + pkin(8) * t293;
t250 = t223 * t304 + (-t201 + t275) * t230 + t258;
t182 = t269 * t243;
t249 = t230 * t182 + t241 * t282 + t252;
t200 = t241 * pkin(4) + pkin(8) * t292;
t248 = t224 * t201 + (-t182 - t200) * t223 + t251;
t247 = t230 * t200 + (-t197 - t304) * t224 + t249;
t174 = -Icges(3,6) * t243 + t241 * t266;
t175 = Icges(3,6) * t241 + t243 * t266;
t176 = -Icges(3,5) * t243 + t241 * t268;
t177 = Icges(3,5) * t241 + t243 * t268;
t208 = Icges(3,2) * t237 + t299;
t209 = Icges(3,1) * t236 + t298;
t244 = (-t175 * t236 + t177 * t237) * V_base(4) + (-t174 * t236 + t176 * t237) * V_base(5) + (-t208 * t236 + t209 * t237) * t230;
t233 = Icges(2,4) * t243;
t222 = rSges(2,1) * t243 - t241 * rSges(2,2);
t220 = t241 * rSges(2,1) + rSges(2,2) * t243;
t216 = Icges(2,1) * t243 - t300;
t215 = Icges(2,1) * t241 + t233;
t214 = -Icges(2,2) * t241 + t233;
t213 = Icges(2,2) * t243 + t300;
t212 = Icges(2,5) * t243 - Icges(2,6) * t241;
t211 = Icges(2,5) * t241 + Icges(2,6) * t243;
t210 = rSges(3,1) * t236 + rSges(3,2) * t237;
t206 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t205 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t204 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t203 = qJD(5) * t228 + t230;
t199 = rSges(4,1) * t228 + rSges(4,2) * t229;
t198 = -rSges(5,2) * t228 - rSges(5,3) * t229;
t184 = t243 * t281 + t224;
t183 = t241 * t281 + t223;
t179 = t241 * rSges(3,3) + t243 * t272;
t178 = -rSges(3,3) * t243 + t241 * t272;
t169 = -rSges(5,1) * t243 + t241 * t270;
t168 = t241 * rSges(5,1) + t243 * t270;
t167 = t241 * rSges(4,3) + t243 * t271;
t166 = -rSges(4,3) * t243 + t241 * t271;
t149 = rSges(6,3) * t228 + (-rSges(6,1) * t240 - rSges(6,2) * t242) * t229;
t147 = V_base(5) * rSges(2,3) - t220 * t230 + t278;
t146 = t222 * t230 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t138 = t220 * V_base(4) - t222 * V_base(5) + V_base(3);
t133 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t293;
t131 = t186 * rSges(6,1) + t185 * rSges(6,2) + rSges(6,3) * t292;
t117 = t210 * V_base(5) + (-t178 - t219) * t230 + t274;
t116 = t230 * t179 + (-pkin(6) - t210) * V_base(4) + t259;
t115 = t178 * V_base(4) + (-t179 - t221) * V_base(5) + t279;
t114 = t199 * t223 + (-t166 + t283) * t230 + t273;
t113 = t230 * t167 - t224 * t199 + t252;
t112 = t166 * t224 - t167 * t223 + t254;
t111 = t198 * t223 + (-t169 + t275) * t230 + t258;
t110 = t230 * t168 + (-t197 - t198) * t224 + t249;
t109 = t169 * t224 + (-t168 - t182) * t223 + t251;
t108 = -t133 * t203 + t149 * t183 + t250;
t107 = t203 * t131 - t184 * t149 + t247;
t106 = -t131 * t183 + t133 * t184 + t248;
t105 = t183 * t284 - t203 * t285 + t243 * t280 + t250;
t104 = -t184 * t284 + t203 * t286 + t241 * t280 + t247;
t103 = qJD(6) * t228 - t183 * t286 + t184 * t285 + t248;
t1 = m(1) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + m(2) * (t138 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + ((t316 * t187 + t315 * t188 + t317 * t293) * t203 + (t321 * t187 + t319 * t188 + t323 * t293) * t184 + (t320 * t187 + t318 * t188 + t322 * t293) * t183) * t183 / 0.2e1 + ((t316 * t185 + t315 * t186 + t317 * t292) * t203 + (t321 * t185 + t319 * t186 + t323 * t292) * t184 + (t320 * t185 + t318 * t186 + t322 * t292) * t183) * t184 / 0.2e1 + (((-t315 * t240 - t316 * t242) * t203 + (-t319 * t240 - t321 * t242) * t184 + (-t318 * t240 - t320 * t242) * t183) * t229 + (t322 * t183 + t323 * t184 + t317 * t203) * t228) * t203 / 0.2e1 + (t314 * t241 - t313 * t243) * t223 / 0.2e1 + (t313 * t241 + t314 * t243) * t224 / 0.2e1 + (t212 * t230 + t253 * t241 + t244 * t243 + (-t241 * t213 + t215 * t243 + Icges(1,4)) * V_base(5) + (-t241 * t214 + t216 * t243 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t211 * t230 + t244 * t241 - t253 * t243 + (t213 * t243 + t241 * t215 + Icges(1,2)) * V_base(5) + (t214 * t243 + t241 * t216 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t174 * t237 + t176 * t236 + t211) * V_base(5) + (t175 * t237 + t177 * t236 + t212) * V_base(4) + (t327 * t228 + t329 * t229) * t224 + (t328 * t228 + t330 * t229) * t223 + (t208 * t237 + t209 * t236 + t325 * t228 - t326 * t229 + Icges(2,3)) * t230) * t230 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
