% Calculate kinetic energy for
% S6RRRPRR5
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:30
% EndTime: 2019-03-09 18:20:35
% DurationCPUTime: 4.46s
% Computational Cost: add. (1858->339), mult. (2033->501), div. (0->0), fcn. (1871->10), ass. (0->175)
t348 = Icges(4,4) + Icges(5,6);
t347 = Icges(4,1) + Icges(5,2);
t346 = -Icges(4,2) - Icges(5,3);
t255 = qJ(2) + qJ(3);
t251 = cos(t255);
t345 = t348 * t251;
t249 = sin(t255);
t344 = t348 * t249;
t343 = Icges(5,4) - Icges(4,5);
t342 = Icges(5,5) - Icges(4,6);
t341 = t249 * t346 + t345;
t340 = t251 * t347 - t344;
t339 = Icges(5,1) + Icges(4,3);
t258 = sin(qJ(1));
t261 = cos(qJ(1));
t338 = t341 * t258 + t342 * t261;
t337 = -t342 * t258 + t341 * t261;
t336 = t340 * t258 + t343 * t261;
t335 = -t343 * t258 + t340 * t261;
t334 = t251 * t346 - t344;
t333 = t249 * t347 + t345;
t332 = t342 * t249 - t343 * t251;
t214 = V_base(5) + (-qJD(2) - qJD(3)) * t261;
t241 = qJD(2) * t258 + V_base(4);
t215 = qJD(3) * t258 + t241;
t245 = V_base(6) + qJD(1);
t331 = (t334 * t249 + t333 * t251) * t245 + (-t337 * t249 + t335 * t251) * t215 + (-t338 * t249 + t336 * t251) * t214;
t330 = (-t343 * t249 - t342 * t251) * t245 + (t339 * t258 + t332 * t261) * t215 + (t332 * t258 - t339 * t261) * t214;
t257 = sin(qJ(2));
t326 = pkin(2) * t257;
t256 = sin(qJ(5));
t325 = pkin(5) * t256;
t324 = pkin(9) * t249;
t260 = cos(qJ(2));
t323 = pkin(2) * t260;
t259 = cos(qJ(5));
t322 = pkin(5) * t259;
t319 = Icges(2,4) * t258;
t318 = Icges(3,4) * t257;
t317 = Icges(3,4) * t260;
t254 = qJ(5) + qJ(6);
t248 = sin(t254);
t312 = t248 * t258;
t311 = t248 * t261;
t250 = cos(t254);
t310 = t250 * t258;
t309 = t250 * t261;
t308 = t251 * t258;
t307 = t251 * t261;
t306 = t256 * t258;
t305 = t256 * t261;
t304 = t258 * t259;
t303 = t259 * t261;
t159 = -pkin(8) * t261 + t258 * t323;
t238 = t258 * pkin(1) - t261 * pkin(7);
t302 = -t159 - t238;
t301 = qJD(4) * t249;
t300 = qJD(5) * t251;
t299 = qJD(6) * t251;
t298 = V_base(5) * pkin(6) + V_base(1);
t290 = pkin(3) * t251 + qJ(4) * t249;
t195 = t290 * t258;
t295 = -t195 + t302;
t240 = -qJD(2) * t261 + V_base(5);
t294 = t240 * t326 + t298;
t221 = qJD(5) * t249 + t245;
t293 = rSges(3,1) * t260 - rSges(3,2) * t257;
t292 = rSges(4,1) * t251 - rSges(4,2) * t249;
t291 = -rSges(5,2) * t251 + rSges(5,3) * t249;
t180 = t261 * t300 + t215;
t289 = Icges(3,1) * t260 - t318;
t287 = -Icges(3,2) * t257 + t317;
t284 = Icges(3,5) * t260 - Icges(3,6) * t257;
t211 = pkin(3) * t249 - qJ(4) * t251;
t280 = t214 * t211 + t261 * t301 + t294;
t239 = t261 * pkin(1) + t258 * pkin(7);
t279 = -V_base(4) * pkin(6) + t245 * t239 + V_base(2);
t278 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t277 = pkin(10) * t251 + t249 * t325;
t179 = t258 * t300 + t214;
t274 = (-Icges(3,3) * t261 + t258 * t284) * t240 + (Icges(3,3) * t258 + t261 * t284) * t241 + (Icges(3,5) * t257 + Icges(3,6) * t260) * t245;
t160 = pkin(8) * t258 + t261 * t323;
t273 = t241 * t159 - t160 * t240 + t278;
t272 = t245 * t160 - t241 * t326 + t279;
t196 = t290 * t261;
t271 = t245 * t196 + t258 * t301 + t272;
t204 = -pkin(4) * t261 + pkin(9) * t308;
t270 = t214 * t324 + (-t204 + t295) * t245 + t280;
t269 = -qJD(4) * t251 + t215 * t195 + t273;
t203 = pkin(4) * t258 + pkin(9) * t307;
t268 = t215 * t204 + (-t196 - t203) * t214 + t269;
t267 = t245 * t203 + (-t211 - t324) * t215 + t271;
t187 = -Icges(3,6) * t261 + t258 * t287;
t188 = Icges(3,6) * t258 + t261 * t287;
t189 = -Icges(3,5) * t261 + t258 * t289;
t190 = Icges(3,5) * t258 + t261 * t289;
t225 = Icges(3,2) * t260 + t318;
t228 = Icges(3,1) * t257 + t317;
t264 = (-t188 * t257 + t190 * t260) * t241 + (-t187 * t257 + t189 * t260) * t240 + (-t225 * t257 + t228 * t260) * t245;
t252 = Icges(2,4) * t261;
t233 = rSges(2,1) * t261 - rSges(2,2) * t258;
t232 = rSges(2,1) * t258 + rSges(2,2) * t261;
t231 = rSges(3,1) * t257 + rSges(3,2) * t260;
t230 = Icges(2,1) * t261 - t319;
t229 = Icges(2,1) * t258 + t252;
t227 = -Icges(2,2) * t258 + t252;
t226 = Icges(2,2) * t261 + t319;
t220 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t219 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t218 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = rSges(4,1) * t249 + rSges(4,2) * t251;
t212 = -rSges(5,2) * t249 - rSges(5,3) * t251;
t201 = t249 * t306 - t303;
t200 = t249 * t304 + t305;
t199 = t249 * t305 + t304;
t198 = t249 * t303 - t306;
t197 = qJD(6) * t249 + t221;
t192 = rSges(3,3) * t258 + t261 * t293;
t191 = -rSges(3,3) * t261 + t258 * t293;
t184 = t249 * t312 - t309;
t183 = t249 * t310 + t311;
t182 = t249 * t311 + t310;
t181 = t249 * t309 - t312;
t178 = -rSges(5,1) * t261 + t258 * t291;
t177 = rSges(5,1) * t258 + t261 * t291;
t176 = rSges(4,3) * t258 + t261 * t292;
t175 = -rSges(4,3) * t261 + t258 * t292;
t174 = pkin(10) * t249 - t251 * t325;
t158 = rSges(6,3) * t249 + (-rSges(6,1) * t256 - rSges(6,2) * t259) * t251;
t157 = Icges(6,5) * t249 + (-Icges(6,1) * t256 - Icges(6,4) * t259) * t251;
t156 = Icges(6,6) * t249 + (-Icges(6,4) * t256 - Icges(6,2) * t259) * t251;
t155 = Icges(6,3) * t249 + (-Icges(6,5) * t256 - Icges(6,6) * t259) * t251;
t154 = V_base(5) * rSges(2,3) - t232 * t245 + t298;
t153 = t233 * t245 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t150 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t149 = rSges(7,3) * t249 + (-rSges(7,1) * t248 - rSges(7,2) * t250) * t251;
t148 = Icges(7,5) * t249 + (-Icges(7,1) * t248 - Icges(7,4) * t250) * t251;
t147 = Icges(7,6) * t249 + (-Icges(7,4) * t248 - Icges(7,2) * t250) * t251;
t146 = Icges(7,3) * t249 + (-Icges(7,5) * t248 - Icges(7,6) * t250) * t251;
t144 = t261 * t299 + t180;
t143 = t258 * t299 + t179;
t140 = t258 * t277 - t261 * t322;
t139 = t258 * t322 + t261 * t277;
t138 = rSges(6,1) * t201 + rSges(6,2) * t200 + rSges(6,3) * t308;
t137 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t307;
t136 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t308;
t135 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t307;
t134 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t308;
t133 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t307;
t132 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t308;
t131 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t307;
t130 = rSges(7,1) * t184 + rSges(7,2) * t183 + rSges(7,3) * t308;
t129 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t307;
t128 = Icges(7,1) * t184 + Icges(7,4) * t183 + Icges(7,5) * t308;
t127 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t307;
t126 = Icges(7,4) * t184 + Icges(7,2) * t183 + Icges(7,6) * t308;
t125 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t307;
t124 = Icges(7,5) * t184 + Icges(7,6) * t183 + Icges(7,3) * t308;
t123 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t307;
t122 = t231 * t240 + (-t191 - t238) * t245 + t298;
t121 = t192 * t245 - t231 * t241 + t279;
t120 = t191 * t241 - t192 * t240 + t278;
t119 = t213 * t214 + (-t175 + t302) * t245 + t294;
t118 = t176 * t245 - t213 * t215 + t272;
t117 = t175 * t215 - t176 * t214 + t273;
t116 = t212 * t214 + (-t178 + t295) * t245 + t280;
t115 = t177 * t245 + (-t211 - t212) * t215 + t271;
t114 = t178 * t215 + (-t177 - t196) * t214 + t269;
t113 = -t138 * t221 + t158 * t179 + t270;
t112 = t137 * t221 - t158 * t180 + t267;
t111 = -t137 * t179 + t138 * t180 + t268;
t110 = -t130 * t197 - t140 * t221 + t143 * t149 + t174 * t179 + t270;
t109 = t129 * t197 + t139 * t221 - t144 * t149 - t174 * t180 + t267;
t108 = -t129 * t143 + t130 * t144 - t139 * t179 + t140 * t180 + t268;
t1 = t180 * ((t131 * t307 + t198 * t133 + t199 * t135) * t180 + (t132 * t307 + t134 * t198 + t136 * t199) * t179 + (t155 * t307 + t156 * t198 + t157 * t199) * t221) / 0.2e1 + t144 * ((t123 * t307 + t181 * t125 + t182 * t127) * t144 + (t124 * t307 + t126 * t181 + t128 * t182) * t143 + (t146 * t307 + t147 * t181 + t148 * t182) * t197) / 0.2e1 + t179 * ((t131 * t308 + t133 * t200 + t135 * t201) * t180 + (t132 * t308 + t200 * t134 + t201 * t136) * t179 + (t155 * t308 + t156 * t200 + t157 * t201) * t221) / 0.2e1 + t143 * ((t123 * t308 + t125 * t183 + t127 * t184) * t144 + (t124 * t308 + t183 * t126 + t184 * t128) * t143 + (t146 * t308 + t147 * t183 + t148 * t184) * t197) / 0.2e1 + t241 * (t274 * t258 + t264 * t261) / 0.2e1 + t240 * (t264 * t258 - t274 * t261) / 0.2e1 + t221 * ((t131 * t180 + t132 * t179 + t155 * t221) * t249 + ((-t133 * t259 - t135 * t256) * t180 + (-t134 * t259 - t136 * t256) * t179 + (-t156 * t259 - t157 * t256) * t221) * t251) / 0.2e1 + t197 * ((t123 * t144 + t124 * t143 + t146 * t197) * t249 + ((-t125 * t250 - t127 * t248) * t144 + (-t126 * t250 - t128 * t248) * t143 + (-t147 * t250 - t148 * t248) * t197) * t251) / 0.2e1 + m(1) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(2) * (t150 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + (t331 * t258 - t330 * t261) * t214 / 0.2e1 + (t330 * t258 + t331 * t261) * t215 / 0.2e1 + ((-t226 * t258 + t229 * t261 + Icges(1,4)) * V_base(5) + (-t227 * t258 + t230 * t261 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t226 * t261 + t229 * t258 + Icges(1,2)) * V_base(5) + (t227 * t261 + t230 * t258 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t188 * t260 + t190 * t257) * t241 + (t187 * t260 + t189 * t257) * t240 + (t335 * t249 + t337 * t251) * t215 + (t336 * t249 + t338 * t251) * t214 + (t225 * t260 + t228 * t257 + t333 * t249 - t334 * t251 + Icges(2,3)) * t245) * t245 / 0.2e1 + V_base(4) * t245 * (Icges(2,5) * t261 - Icges(2,6) * t258) + V_base(5) * t245 * (Icges(2,5) * t258 + Icges(2,6) * t261) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
