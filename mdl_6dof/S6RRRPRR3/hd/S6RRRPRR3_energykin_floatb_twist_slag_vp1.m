% Calculate kinetic energy for
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:48
% EndTime: 2019-03-09 18:10:52
% DurationCPUTime: 4.18s
% Computational Cost: add. (1977->330), mult. (2393->479), div. (0->0), fcn. (2409->10), ass. (0->164)
t333 = Icges(4,4) - Icges(5,5);
t332 = Icges(4,1) + Icges(5,1);
t331 = Icges(4,2) + Icges(5,3);
t251 = qJ(2) + qJ(3);
t248 = cos(t251);
t330 = t333 * t248;
t247 = sin(t251);
t329 = t333 * t247;
t328 = Icges(5,4) + Icges(4,5);
t327 = Icges(4,6) - Icges(5,6);
t326 = t331 * t247 - t330;
t325 = t332 * t248 - t329;
t324 = Icges(5,2) + Icges(4,3);
t255 = sin(qJ(1));
t258 = cos(qJ(1));
t323 = t326 * t255 + t327 * t258;
t322 = -t327 * t255 + t326 * t258;
t321 = t325 * t255 - t328 * t258;
t320 = t328 * t255 + t325 * t258;
t319 = -t331 * t248 - t329;
t318 = t332 * t247 + t330;
t317 = -t327 * t247 + t328 * t248;
t215 = V_base(5) + (-qJD(2) - qJD(3)) * t258;
t241 = qJD(2) * t255 + V_base(4);
t216 = qJD(3) * t255 + t241;
t243 = V_base(6) + qJD(1);
t316 = (t247 * t319 + t248 * t318) * t243 + (t247 * t322 + t248 * t320) * t216 + (t247 * t323 + t248 * t321) * t215;
t315 = (t328 * t247 + t327 * t248) * t243 + (t255 * t324 + t317 * t258) * t216 + (t317 * t255 - t258 * t324) * t215;
t311 = cos(qJ(5));
t254 = sin(qJ(2));
t310 = pkin(2) * t254;
t309 = pkin(4) * t247;
t257 = cos(qJ(2));
t308 = pkin(2) * t257;
t306 = Icges(2,4) * t255;
t305 = Icges(3,4) * t254;
t304 = Icges(3,4) * t257;
t299 = t248 * t255;
t298 = t248 * t258;
t163 = -pkin(8) * t258 + t255 * t308;
t238 = t255 * pkin(1) - t258 * pkin(7);
t297 = -t163 - t238;
t296 = qJD(4) * t247;
t295 = V_base(5) * pkin(6) + V_base(1);
t285 = pkin(3) * t248 + qJ(4) * t247;
t197 = t285 * t255;
t292 = -t197 + t297;
t291 = t247 * t311;
t240 = -qJD(2) * t258 + V_base(5);
t290 = t240 * t310 + t295;
t204 = pkin(4) * t299 + pkin(9) * t258;
t289 = -t204 + t292;
t288 = rSges(3,1) * t257 - rSges(3,2) * t254;
t287 = rSges(4,1) * t248 - rSges(4,2) * t247;
t286 = rSges(5,1) * t248 + rSges(5,3) * t247;
t284 = Icges(3,1) * t257 - t305;
t281 = -Icges(3,2) * t254 + t304;
t278 = Icges(3,5) * t257 - Icges(3,6) * t254;
t212 = pkin(3) * t247 - qJ(4) * t248;
t275 = t215 * t212 + t258 * t296 + t290;
t239 = t258 * pkin(1) + t255 * pkin(7);
t274 = -V_base(4) * pkin(6) + t243 * t239 + V_base(2);
t253 = sin(qJ(5));
t201 = t247 * t253 + t248 * t311;
t273 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t272 = t215 * t309 + t275;
t200 = -qJD(5) * t255 + t216;
t199 = qJD(5) * t258 + t215;
t269 = (-Icges(3,3) * t258 + t255 * t278) * t240 + (Icges(3,3) * t255 + t258 * t278) * t241 + (Icges(3,5) * t254 + Icges(3,6) * t257) * t243;
t164 = pkin(8) * t255 + t258 * t308;
t268 = t241 * t163 - t164 * t240 + t273;
t267 = t243 * t164 - t241 * t310 + t274;
t198 = t285 * t258;
t266 = t243 * t198 + t255 * t296 + t267;
t265 = -qJD(4) * t248 + t216 * t197 + t268;
t205 = pkin(4) * t298 - pkin(9) * t255;
t264 = t216 * t204 + (-t198 - t205) * t215 + t265;
t263 = t243 * t205 + (-t212 - t309) * t216 + t266;
t189 = -Icges(3,6) * t258 + t255 * t281;
t190 = Icges(3,6) * t255 + t258 * t281;
t191 = -Icges(3,5) * t258 + t255 * t284;
t192 = Icges(3,5) * t255 + t258 * t284;
t227 = Icges(3,2) * t257 + t305;
t230 = Icges(3,1) * t254 + t304;
t260 = (-t190 * t254 + t192 * t257) * t241 + (-t189 * t254 + t191 * t257) * t240 + (-t227 * t254 + t230 * t257) * t243;
t256 = cos(qJ(6));
t252 = sin(qJ(6));
t249 = Icges(2,4) * t258;
t235 = rSges(2,1) * t258 - rSges(2,2) * t255;
t234 = rSges(2,1) * t255 + rSges(2,2) * t258;
t233 = rSges(3,1) * t254 + rSges(3,2) * t257;
t232 = Icges(2,1) * t258 - t306;
t231 = Icges(2,1) * t255 + t249;
t229 = -Icges(2,2) * t255 + t249;
t228 = Icges(2,2) * t258 + t306;
t221 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t220 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t219 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t214 = rSges(4,1) * t247 + rSges(4,2) * t248;
t213 = rSges(5,1) * t247 - rSges(5,3) * t248;
t202 = -t248 * t253 + t291;
t194 = rSges(3,3) * t255 + t258 * t288;
t193 = -rSges(3,3) * t258 + t255 * t288;
t186 = t201 * t258;
t185 = t253 * t298 - t258 * t291;
t184 = t201 * t255;
t183 = t253 * t299 - t255 * t291;
t182 = rSges(4,3) * t255 + t258 * t287;
t181 = rSges(5,2) * t255 + t258 * t286;
t180 = -rSges(4,3) * t258 + t255 * t287;
t179 = -rSges(5,2) * t258 + t255 * t286;
t178 = qJD(6) * t201 + t243;
t162 = V_base(5) * rSges(2,3) - t234 * t243 + t295;
t161 = t235 * t243 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t158 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t156 = t186 * t256 - t252 * t255;
t155 = -t186 * t252 - t255 * t256;
t154 = t184 * t256 + t252 * t258;
t153 = -t184 * t252 + t256 * t258;
t150 = pkin(5) * t202 + pkin(10) * t201;
t149 = rSges(6,1) * t202 - rSges(6,2) * t201;
t148 = Icges(6,1) * t202 - Icges(6,4) * t201;
t147 = Icges(6,4) * t202 - Icges(6,2) * t201;
t146 = Icges(6,5) * t202 - Icges(6,6) * t201;
t145 = qJD(6) * t185 + t200;
t144 = qJD(6) * t183 + t199;
t143 = pkin(5) * t186 + pkin(10) * t185;
t142 = pkin(5) * t184 + pkin(10) * t183;
t141 = rSges(6,1) * t186 - rSges(6,2) * t185 - rSges(6,3) * t255;
t140 = rSges(6,1) * t184 - rSges(6,2) * t183 + rSges(6,3) * t258;
t139 = Icges(6,1) * t186 - Icges(6,4) * t185 - Icges(6,5) * t255;
t138 = Icges(6,1) * t184 - Icges(6,4) * t183 + Icges(6,5) * t258;
t137 = Icges(6,4) * t186 - Icges(6,2) * t185 - Icges(6,6) * t255;
t136 = Icges(6,4) * t184 - Icges(6,2) * t183 + Icges(6,6) * t258;
t135 = Icges(6,5) * t186 - Icges(6,6) * t185 - Icges(6,3) * t255;
t134 = Icges(6,5) * t184 - Icges(6,6) * t183 + Icges(6,3) * t258;
t133 = rSges(7,3) * t201 + (rSges(7,1) * t256 - rSges(7,2) * t252) * t202;
t132 = Icges(7,5) * t201 + (Icges(7,1) * t256 - Icges(7,4) * t252) * t202;
t131 = Icges(7,6) * t201 + (Icges(7,4) * t256 - Icges(7,2) * t252) * t202;
t130 = Icges(7,3) * t201 + (Icges(7,5) * t256 - Icges(7,6) * t252) * t202;
t129 = t233 * t240 + (-t193 - t238) * t243 + t295;
t128 = t194 * t243 - t233 * t241 + t274;
t127 = t193 * t241 - t194 * t240 + t273;
t126 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t185;
t125 = rSges(7,1) * t154 + rSges(7,2) * t153 + rSges(7,3) * t183;
t124 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t185;
t123 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t183;
t122 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t185;
t121 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t183;
t120 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t185;
t119 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t183;
t118 = t214 * t215 + (-t180 + t297) * t243 + t290;
t117 = t182 * t243 - t214 * t216 + t267;
t116 = t180 * t216 - t182 * t215 + t268;
t115 = t213 * t215 + (-t179 + t292) * t243 + t275;
t114 = t181 * t243 + (-t212 - t213) * t216 + t266;
t113 = t179 * t216 + (-t181 - t198) * t215 + t265;
t112 = t149 * t199 + (-t140 + t289) * t243 + t272;
t111 = t141 * t243 - t149 * t200 + t263;
t110 = t140 * t200 - t141 * t199 + t264;
t109 = -t125 * t178 + t133 * t144 + t150 * t199 + (-t142 + t289) * t243 + t272;
t108 = t126 * t178 - t133 * t145 + t143 * t243 - t150 * t200 + t263;
t107 = t125 * t145 - t126 * t144 + t142 * t200 - t143 * t199 + t264;
t1 = t241 * (t269 * t255 + t260 * t258) / 0.2e1 + t240 * (t260 * t255 - t269 * t258) / 0.2e1 + t178 * ((t119 * t144 + t120 * t145 + t130 * t178) * t201 + ((-t122 * t252 + t124 * t256) * t145 + (-t121 * t252 + t123 * t256) * t144 + (-t131 * t252 + t132 * t256) * t178) * t202) / 0.2e1 + t199 * ((t135 * t258 - t137 * t183 + t139 * t184) * t200 + (t258 * t134 - t183 * t136 + t184 * t138) * t199 + (t146 * t258 - t147 * t183 + t148 * t184) * t243) / 0.2e1 + t200 * ((-t255 * t135 - t185 * t137 + t186 * t139) * t200 + (-t134 * t255 - t136 * t185 + t138 * t186) * t199 + (-t146 * t255 - t147 * t185 + t148 * t186) * t243) / 0.2e1 + m(1) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + t144 * ((t120 * t183 + t122 * t153 + t124 * t154) * t145 + (t183 * t119 + t153 * t121 + t154 * t123) * t144 + (t130 * t183 + t131 * t153 + t132 * t154) * t178) / 0.2e1 + t145 * ((t185 * t120 + t155 * t122 + t156 * t124) * t145 + (t119 * t185 + t121 * t155 + t123 * t156) * t144 + (t130 * t185 + t131 * t155 + t132 * t156) * t178) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(2) * (t158 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + (t316 * t255 - t315 * t258) * t215 / 0.2e1 + (t315 * t255 + t316 * t258) * t216 / 0.2e1 + ((-t228 * t255 + t231 * t258 + Icges(1,4)) * V_base(5) + (-t229 * t255 + t232 * t258 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t228 * t258 + t231 * t255 + Icges(1,2)) * V_base(5) + (t229 * t258 + t232 * t255 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t190 * t257 + t192 * t254) * t241 + (t189 * t257 + t191 * t254) * t240 + (-t137 * t201 + t139 * t202) * t200 + (-t136 * t201 + t138 * t202) * t199 + (t247 * t320 - t248 * t322) * t216 + (t247 * t321 - t248 * t323) * t215 + (-t147 * t201 + t148 * t202 + t227 * t257 + t230 * t254 + t318 * t247 - t319 * t248 + Icges(2,3)) * t243) * t243 / 0.2e1 + t243 * V_base(4) * (Icges(2,5) * t258 - Icges(2,6) * t255) + V_base(5) * t243 * (Icges(2,5) * t255 + Icges(2,6) * t258) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
