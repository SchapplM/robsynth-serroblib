% Calculate kinetic energy for
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:34
% EndTime: 2019-03-09 08:29:37
% DurationCPUTime: 3.52s
% Computational Cost: add. (1624->288), mult. (1944->389), div. (0->0), fcn. (1782->8), ass. (0->152)
t347 = Icges(4,4) + Icges(5,6);
t346 = Icges(4,1) + Icges(5,2);
t345 = -Icges(4,2) - Icges(5,3);
t235 = qJ(2) + pkin(9);
t229 = cos(t235);
t344 = t347 * t229;
t228 = sin(t235);
t343 = t347 * t228;
t342 = Icges(5,4) - Icges(4,5);
t341 = Icges(5,5) - Icges(4,6);
t340 = t228 * t345 + t344;
t339 = t229 * t346 - t343;
t338 = Icges(6,1) + Icges(7,1);
t337 = Icges(6,4) + Icges(7,4);
t336 = Icges(7,5) + Icges(6,5);
t335 = Icges(6,2) + Icges(7,2);
t334 = Icges(7,6) + Icges(6,6);
t333 = Icges(7,3) + Icges(6,3);
t240 = sin(qJ(1));
t243 = cos(qJ(1));
t332 = t240 * t340 + t243 * t341;
t331 = -t240 * t341 + t243 * t340;
t330 = t240 * t339 + t243 * t342;
t329 = -t240 * t342 + t243 * t339;
t328 = t229 * t345 - t343;
t327 = t228 * t346 + t344;
t326 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t239 = sin(qJ(2));
t242 = cos(qJ(2));
t325 = Icges(3,5) * t242 - Icges(3,6) * t239 + t228 * t341 - t229 * t342;
t241 = cos(qJ(5));
t288 = t241 * t243;
t238 = sin(qJ(5));
t290 = t240 * t238;
t185 = t228 * t288 - t290;
t289 = t240 * t241;
t291 = t238 * t243;
t186 = t228 * t291 + t289;
t292 = t229 * t243;
t324 = t185 * t334 + t186 * t336 + t292 * t333;
t187 = t228 * t289 + t291;
t188 = t228 * t290 - t288;
t293 = t229 * t240;
t323 = t187 * t334 + t188 * t336 + t293 * t333;
t322 = t185 * t335 + t186 * t337 + t292 * t334;
t321 = t187 * t335 + t188 * t337 + t293 * t334;
t320 = t185 * t337 + t186 * t338 + t292 * t336;
t319 = t187 * t337 + t188 * t338 + t293 * t336;
t318 = (-t238 * t336 - t241 * t334) * t229 + t333 * t228;
t317 = (-t238 * t337 - t241 * t335) * t229 + t334 * t228;
t316 = (-t238 * t338 - t241 * t337) * t229 + t336 * t228;
t298 = Icges(3,4) * t242;
t267 = -Icges(3,2) * t239 + t298;
t174 = -Icges(3,6) * t243 + t240 * t267;
t175 = Icges(3,6) * t240 + t243 * t267;
t299 = Icges(3,4) * t239;
t269 = Icges(3,1) * t242 - t299;
t176 = -Icges(3,5) * t243 + t240 * t269;
t177 = Icges(3,5) * t240 + t243 * t269;
t211 = Icges(3,2) * t242 + t299;
t214 = Icges(3,1) * t239 + t298;
t224 = -qJD(2) * t243 + V_base(5);
t225 = qJD(2) * t240 + V_base(4);
t230 = V_base(6) + qJD(1);
t315 = (-t211 * t239 + t214 * t242 + t228 * t328 + t229 * t327) * t230 + (-t175 * t239 + t177 * t242 - t228 * t331 + t229 * t329) * t225 + (-t174 * t239 + t176 * t242 - t228 * t332 + t229 * t330) * t224;
t314 = (Icges(3,5) * t239 + Icges(3,6) * t242 - t228 * t342 - t229 * t341) * t230 + (t240 * t326 + t243 * t325) * t225 + (t240 * t325 - t243 * t326) * t224;
t307 = pkin(2) * t239;
t306 = pkin(5) * t238;
t305 = pkin(8) * t228;
t304 = pkin(2) * t242;
t303 = pkin(5) * t241;
t300 = Icges(2,4) * t240;
t257 = qJ(6) * t229 + t228 * t306;
t287 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t292 + t240 * t303 + t243 * t257;
t286 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t293 + t240 * t257 - t243 * t303;
t285 = (-rSges(7,1) * t238 - rSges(7,2) * t241 - t306) * t229 + (rSges(7,3) + qJ(6)) * t228;
t152 = -qJ(3) * t243 + t240 * t304;
t222 = pkin(1) * t240 - pkin(7) * t243;
t284 = -t152 - t222;
t153 = qJ(3) * t240 + t243 * t304;
t270 = pkin(3) * t229 + qJ(4) * t228;
t180 = t270 * t243;
t283 = -t153 - t180;
t282 = qJD(4) * t228;
t281 = qJD(5) * t229;
t280 = qJD(6) * t229;
t279 = V_base(5) * pkin(6) + V_base(1);
t179 = t270 * t240;
t276 = -t179 + t284;
t197 = pkin(3) * t228 - qJ(4) * t229;
t275 = -t197 - t307;
t274 = qJD(3) * t240 + t224 * t307 + t279;
t273 = rSges(3,1) * t242 - rSges(3,2) * t239;
t272 = rSges(4,1) * t229 - rSges(4,2) * t228;
t271 = -rSges(5,2) * t229 + rSges(5,3) * t228;
t223 = pkin(1) * t243 + pkin(7) * t240;
t260 = -V_base(4) * pkin(6) + t223 * t230 + V_base(2);
t259 = t222 * V_base(4) - t223 * V_base(5) + V_base(3);
t258 = t197 * t224 + t243 * t282 + t274;
t256 = t152 * t225 + t259;
t252 = -qJD(3) * t243 + t153 * t230 + t260;
t251 = -qJD(4) * t229 + t179 * t225 + t256;
t250 = t180 * t230 + t240 * t282 + t252;
t201 = -pkin(4) * t243 + pkin(8) * t293;
t249 = t224 * t305 + (-t201 + t276) * t230 + t258;
t200 = pkin(4) * t240 + pkin(8) * t292;
t248 = t225 * t201 + (-t200 + t283) * t224 + t251;
t247 = t230 * t200 + (t275 - t305) * t225 + t250;
t233 = Icges(2,4) * t243;
t221 = rSges(2,1) * t243 - rSges(2,2) * t240;
t220 = rSges(2,1) * t240 + rSges(2,2) * t243;
t219 = rSges(3,1) * t239 + rSges(3,2) * t242;
t216 = Icges(2,1) * t243 - t300;
t215 = Icges(2,1) * t240 + t233;
t213 = -Icges(2,2) * t240 + t233;
t212 = Icges(2,2) * t243 + t300;
t207 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t206 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t205 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t204 = qJD(5) * t228 + t230;
t199 = rSges(4,1) * t228 + rSges(4,2) * t229;
t198 = -rSges(5,2) * t228 - rSges(5,3) * t229;
t184 = t243 * t281 + t225;
t183 = t240 * t281 + t224;
t182 = rSges(3,3) * t240 + t243 * t273;
t181 = -rSges(3,3) * t243 + t240 * t273;
t169 = -rSges(5,1) * t243 + t240 * t271;
t168 = rSges(5,1) * t240 + t243 * t271;
t167 = rSges(4,3) * t240 + t243 * t272;
t166 = -rSges(4,3) * t243 + t240 * t272;
t149 = rSges(6,3) * t228 + (-rSges(6,1) * t238 - rSges(6,2) * t241) * t229;
t147 = V_base(5) * rSges(2,3) - t220 * t230 + t279;
t146 = t221 * t230 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t138 = t220 * V_base(4) - t221 * V_base(5) + V_base(3);
t133 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t293;
t131 = rSges(6,1) * t186 + rSges(6,2) * t185 + rSges(6,3) * t292;
t117 = t219 * t224 + (-t181 - t222) * t230 + t279;
t116 = t182 * t230 - t219 * t225 + t260;
t115 = t181 * t225 - t182 * t224 + t259;
t114 = t199 * t224 + (-t166 + t284) * t230 + t274;
t113 = t230 * t167 + (-t199 - t307) * t225 + t252;
t112 = t166 * t225 + (-t153 - t167) * t224 + t256;
t111 = t198 * t224 + (-t169 + t276) * t230 + t258;
t110 = t230 * t168 + (-t198 + t275) * t225 + t250;
t109 = t169 * t225 + (-t168 + t283) * t224 + t251;
t108 = -t133 * t204 + t149 * t183 + t249;
t107 = t131 * t204 - t149 * t184 + t247;
t106 = -t131 * t183 + t133 * t184 + t248;
t105 = t183 * t285 - t204 * t286 + t243 * t280 + t249;
t104 = -t184 * t285 + t204 * t287 + t240 * t280 + t247;
t103 = qJD(6) * t228 - t183 * t287 + t184 * t286 + t248;
t1 = m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(2) * (t138 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + ((t187 * t317 + t188 * t316 + t293 * t318) * t204 + (t187 * t322 + t188 * t320 + t293 * t324) * t184 + (t321 * t187 + t319 * t188 + t323 * t293) * t183) * t183 / 0.2e1 + ((t185 * t317 + t186 * t316 + t292 * t318) * t204 + (t322 * t185 + t320 * t186 + t324 * t292) * t184 + (t185 * t321 + t186 * t319 + t292 * t323) * t183) * t184 / 0.2e1 + (((-t238 * t316 - t241 * t317) * t204 + (-t238 * t320 - t241 * t322) * t184 + (-t238 * t319 - t241 * t321) * t183) * t229 + (t183 * t323 + t184 * t324 + t204 * t318) * t228) * t204 / 0.2e1 + ((-t212 * t240 + t215 * t243 + Icges(1,4)) * V_base(5) + (-t240 * t213 + t216 * t243 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t212 * t243 + t240 * t215 + Icges(1,2)) * V_base(5) + (t213 * t243 + t216 * t240 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t315 * t240 - t314 * t243) * t224 / 0.2e1 + (t314 * t240 + t315 * t243) * t225 / 0.2e1 + ((t175 * t242 + t177 * t239 + t228 * t329 + t229 * t331) * t225 + (t174 * t242 + t176 * t239 + t228 * t330 + t229 * t332) * t224 + (t211 * t242 + t214 * t239 + t228 * t327 - t328 * t229 + Icges(2,3)) * t230) * t230 / 0.2e1 + t230 * V_base(4) * (Icges(2,5) * t243 - Icges(2,6) * t240) + V_base(5) * t230 * (Icges(2,5) * t240 + Icges(2,6) * t243) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
