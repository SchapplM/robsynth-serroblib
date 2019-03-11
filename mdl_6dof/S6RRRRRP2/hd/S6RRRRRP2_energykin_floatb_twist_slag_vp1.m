% Calculate kinetic energy for
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:36
% EndTime: 2019-03-10 01:00:39
% DurationCPUTime: 2.73s
% Computational Cost: add. (2197->309), mult. (2080->454), div. (0->0), fcn. (1928->10), ass. (0->166)
t334 = Icges(6,1) + Icges(7,1);
t333 = -Icges(6,4) + Icges(7,5);
t332 = Icges(7,4) + Icges(6,5);
t331 = Icges(6,2) + Icges(7,3);
t330 = -Icges(7,6) + Icges(6,6);
t329 = -Icges(6,3) - Icges(7,2);
t328 = rSges(7,1) + pkin(5);
t327 = rSges(7,3) + qJ(6);
t247 = qJ(2) + qJ(3);
t243 = qJ(4) + t247;
t234 = cos(t243);
t251 = cos(qJ(5));
t253 = cos(qJ(1));
t297 = t251 * t253;
t248 = sin(qJ(5));
t250 = sin(qJ(1));
t300 = t248 * t250;
t191 = t234 * t300 + t297;
t298 = t250 * t251;
t299 = t248 * t253;
t192 = t234 * t298 - t299;
t233 = sin(t243);
t302 = t233 * t250;
t326 = t331 * t191 + t333 * t192 - t330 * t302;
t193 = t234 * t299 - t298;
t194 = t234 * t297 + t300;
t301 = t233 * t253;
t325 = t331 * t193 + t333 * t194 - t330 * t301;
t324 = -t330 * t191 + t332 * t192 - t329 * t302;
t323 = -t330 * t193 + t332 * t194 - t329 * t301;
t322 = t333 * t191 + t334 * t192 + t332 * t302;
t321 = t333 * t193 + t334 * t194 + t332 * t301;
t320 = t330 * t234 + (t331 * t248 + t333 * t251) * t233;
t319 = t329 * t234 + (-t330 * t248 + t332 * t251) * t233;
t318 = -t332 * t234 + (t333 * t248 + t334 * t251) * t233;
t249 = sin(qJ(2));
t313 = pkin(2) * t249;
t240 = sin(t247);
t312 = pkin(3) * t240;
t252 = cos(qJ(2));
t311 = t252 * pkin(2);
t309 = Icges(2,4) * t250;
t308 = Icges(3,4) * t249;
t307 = Icges(3,4) * t252;
t306 = Icges(4,4) * t240;
t241 = cos(t247);
t305 = Icges(4,4) * t241;
t304 = Icges(5,4) * t233;
t303 = Icges(5,4) * t234;
t296 = rSges(7,2) * t302 + t327 * t191 + t328 * t192;
t295 = rSges(7,2) * t301 + t327 * t193 + t328 * t194;
t294 = -rSges(7,2) * t234 + (t327 * t248 + t328 * t251) * t233;
t169 = -pkin(8) * t253 + t250 * t311;
t228 = t250 * pkin(1) - t253 * pkin(7);
t293 = -t169 - t228;
t292 = pkin(3) * t241;
t290 = qJD(5) * t233;
t289 = -qJD(2) - qJD(3);
t288 = V_base(5) * pkin(6) + V_base(1);
t142 = -pkin(9) * t253 + t250 * t292;
t285 = -t142 + t293;
t231 = qJD(2) * t250 + V_base(4);
t236 = V_base(6) + qJD(1);
t230 = -qJD(2) * t253 + V_base(5);
t284 = t230 * t313 + t288;
t208 = qJD(3) * t250 + t231;
t207 = t253 * t289 + V_base(5);
t283 = t207 * t312 + t284;
t282 = pkin(4) * t234 + pkin(10) * t233;
t281 = rSges(3,1) * t252 - rSges(3,2) * t249;
t280 = rSges(4,1) * t241 - rSges(4,2) * t240;
t279 = rSges(5,1) * t234 - rSges(5,2) * t233;
t196 = qJD(4) * t250 + t208;
t278 = Icges(3,1) * t252 - t308;
t277 = Icges(4,1) * t241 - t306;
t276 = Icges(5,1) * t234 - t304;
t275 = -Icges(3,2) * t249 + t307;
t274 = -Icges(4,2) * t240 + t305;
t273 = -Icges(5,2) * t233 + t303;
t272 = Icges(3,5) * t252 - Icges(3,6) * t249;
t271 = Icges(4,5) * t241 - Icges(4,6) * t240;
t270 = Icges(5,5) * t234 - Icges(5,6) * t233;
t229 = t253 * pkin(1) + t250 * pkin(7);
t269 = -V_base(4) * pkin(6) + t236 * t229 + V_base(2);
t268 = V_base(4) * t228 - t229 * V_base(5) + V_base(3);
t195 = V_base(5) + (-qJD(4) + t289) * t253;
t267 = (-Icges(5,3) * t253 + t250 * t270) * t195 + (Icges(5,3) * t250 + t253 * t270) * t196 + (Icges(5,5) * t233 + Icges(5,6) * t234) * t236;
t266 = (-Icges(4,3) * t253 + t250 * t271) * t207 + (Icges(4,3) * t250 + t253 * t271) * t208 + (Icges(4,5) * t240 + Icges(4,6) * t241) * t236;
t265 = (-Icges(3,3) * t253 + t250 * t272) * t230 + (Icges(3,3) * t250 + t253 * t272) * t231 + (Icges(3,5) * t249 + Icges(3,6) * t252) * t236;
t170 = pkin(8) * t250 + t253 * t311;
t264 = t231 * t169 - t170 * t230 + t268;
t263 = t236 * t170 - t231 * t313 + t269;
t180 = t282 * t250;
t201 = pkin(4) * t233 - pkin(10) * t234;
t262 = t195 * t201 + (-t180 + t285) * t236 + t283;
t143 = pkin(9) * t250 + t253 * t292;
t261 = t208 * t142 - t143 * t207 + t264;
t260 = t236 * t143 - t208 * t312 + t263;
t181 = t282 * t253;
t259 = t196 * t180 - t181 * t195 + t261;
t258 = t236 * t181 - t196 * t201 + t260;
t163 = -Icges(5,6) * t253 + t250 * t273;
t164 = Icges(5,6) * t250 + t253 * t273;
t165 = -Icges(5,5) * t253 + t250 * t276;
t166 = Icges(5,5) * t250 + t253 * t276;
t198 = Icges(5,2) * t234 + t304;
t199 = Icges(5,1) * t233 + t303;
t257 = (-t164 * t233 + t166 * t234) * t196 + (-t163 * t233 + t165 * t234) * t195 + (-t198 * t233 + t199 * t234) * t236;
t173 = -Icges(4,6) * t253 + t250 * t274;
t174 = Icges(4,6) * t250 + t253 * t274;
t175 = -Icges(4,5) * t253 + t250 * t277;
t176 = Icges(4,5) * t250 + t253 * t277;
t204 = Icges(4,2) * t241 + t306;
t205 = Icges(4,1) * t240 + t305;
t256 = (-t174 * t240 + t176 * t241) * t208 + (-t173 * t240 + t175 * t241) * t207 + (-t204 * t240 + t205 * t241) * t236;
t184 = -Icges(3,6) * t253 + t250 * t275;
t185 = Icges(3,6) * t250 + t253 * t275;
t186 = -Icges(3,5) * t253 + t250 * t278;
t187 = Icges(3,5) * t250 + t253 * t278;
t219 = Icges(3,2) * t252 + t308;
t222 = Icges(3,1) * t249 + t307;
t255 = (-t185 * t249 + t187 * t252) * t231 + (-t184 * t249 + t186 * t252) * t230 + (-t219 * t249 + t222 * t252) * t236;
t242 = Icges(2,4) * t253;
t227 = rSges(2,1) * t253 - rSges(2,2) * t250;
t226 = rSges(2,1) * t250 + rSges(2,2) * t253;
t225 = rSges(3,1) * t249 + rSges(3,2) * t252;
t224 = Icges(2,1) * t253 - t309;
t223 = Icges(2,1) * t250 + t242;
t221 = -Icges(2,2) * t250 + t242;
t220 = Icges(2,2) * t253 + t309;
t215 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t214 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t213 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t209 = -qJD(5) * t234 + t236;
t206 = rSges(4,1) * t240 + rSges(4,2) * t241;
t200 = rSges(5,1) * t233 + rSges(5,2) * t234;
t189 = rSges(3,3) * t250 + t253 * t281;
t188 = -rSges(3,3) * t253 + t250 * t281;
t178 = rSges(4,3) * t250 + t253 * t280;
t177 = -rSges(4,3) * t253 + t250 * t280;
t168 = rSges(5,3) * t250 + t253 * t279;
t167 = -rSges(5,3) * t253 + t250 * t279;
t159 = V_base(5) * rSges(2,3) - t226 * t236 + t288;
t158 = t227 * t236 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t157 = t253 * t290 + t196;
t156 = t250 * t290 + t195;
t155 = t226 * V_base(4) - t227 * V_base(5) + V_base(3);
t154 = -rSges(6,3) * t234 + (rSges(6,1) * t251 - rSges(6,2) * t248) * t233;
t136 = rSges(6,1) * t194 - rSges(6,2) * t193 + rSges(6,3) * t301;
t134 = rSges(6,1) * t192 - rSges(6,2) * t191 + rSges(6,3) * t302;
t120 = t225 * t230 + (-t188 - t228) * t236 + t288;
t119 = t189 * t236 - t225 * t231 + t269;
t118 = t188 * t231 - t189 * t230 + t268;
t117 = t206 * t207 + (-t177 + t293) * t236 + t284;
t116 = t178 * t236 - t206 * t208 + t263;
t115 = t177 * t208 - t178 * t207 + t264;
t114 = t195 * t200 + (-t167 + t285) * t236 + t283;
t113 = t168 * t236 - t196 * t200 + t260;
t112 = t167 * t196 - t168 * t195 + t261;
t111 = -t134 * t209 + t154 * t156 + t262;
t110 = t136 * t209 - t154 * t157 + t258;
t109 = qJD(6) * t193 + t156 * t294 - t209 * t296 + t262;
t108 = qJD(6) * t191 - t157 * t294 + t209 * t295 + t258;
t107 = t134 * t157 - t136 * t156 + t259;
t106 = qJD(6) * t233 * t248 - t156 * t295 + t157 * t296 + t259;
t1 = m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(3) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t231 * (t265 * t250 + t255 * t253) / 0.2e1 + t230 * (t255 * t250 - t265 * t253) / 0.2e1 + t208 * (t266 * t250 + t256 * t253) / 0.2e1 + t207 * (t256 * t250 - t266 * t253) / 0.2e1 + t196 * (t267 * t250 + t257 * t253) / 0.2e1 + t195 * (t257 * t250 - t267 * t253) / 0.2e1 + ((t191 * t320 + t192 * t318 + t302 * t319) * t209 + (t191 * t325 + t192 * t321 + t302 * t323) * t157 + (t326 * t191 + t322 * t192 + t324 * t302) * t156) * t156 / 0.2e1 + ((t193 * t320 + t194 * t318 + t301 * t319) * t209 + (t325 * t193 + t321 * t194 + t323 * t301) * t157 + (t193 * t326 + t322 * t194 + t324 * t301) * t156) * t157 / 0.2e1 + ((-t156 * t324 - t157 * t323 - t209 * t319) * t234 + ((t248 * t320 + t251 * t318) * t209 + (t248 * t325 + t251 * t321) * t157 + (t248 * t326 + t322 * t251) * t156) * t233) * t209 / 0.2e1 + ((-t220 * t250 + t223 * t253 + Icges(1,4)) * V_base(5) + (-t221 * t250 + t224 * t253 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t220 * t253 + t223 * t250 + Icges(1,2)) * V_base(5) + (t221 * t253 + t224 * t250 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t185 * t252 + t187 * t249) * t231 + (t184 * t252 + t186 * t249) * t230 + (t174 * t241 + t176 * t240) * t208 + (t173 * t241 + t175 * t240) * t207 + (t164 * t234 + t166 * t233) * t196 + (t163 * t234 + t165 * t233) * t195 + (t198 * t234 + t199 * t233 + t204 * t241 + t205 * t240 + t219 * t252 + t222 * t249 + Icges(2,3)) * t236) * t236 / 0.2e1 + t236 * V_base(4) * (Icges(2,5) * t253 - Icges(2,6) * t250) + t236 * V_base(5) * (Icges(2,5) * t250 + Icges(2,6) * t253) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
