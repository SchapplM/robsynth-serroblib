% Calculate kinetic energy for
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:38:55
% EndTime: 2019-03-09 11:38:57
% DurationCPUTime: 2.84s
% Computational Cost: add. (2175->311), mult. (2025->439), div. (0->0), fcn. (1863->10), ass. (0->162)
t338 = Icges(6,1) + Icges(7,1);
t337 = Icges(6,4) + Icges(7,4);
t336 = -Icges(7,5) - Icges(6,5);
t335 = Icges(6,2) + Icges(7,2);
t334 = -Icges(7,6) - Icges(6,6);
t333 = Icges(3,3) + Icges(4,3);
t332 = -Icges(7,3) - Icges(6,3);
t243 = qJ(2) + pkin(10);
t232 = sin(t243);
t233 = cos(t243);
t247 = sin(qJ(2));
t250 = cos(qJ(2));
t331 = Icges(3,5) * t250 + Icges(4,5) * t233 - Icges(3,6) * t247 - Icges(4,6) * t232;
t234 = qJ(4) + t243;
t229 = cos(t234);
t249 = cos(qJ(5));
t251 = cos(qJ(1));
t295 = t249 * t251;
t246 = sin(qJ(5));
t248 = sin(qJ(1));
t297 = t248 * t246;
t188 = -t229 * t297 - t295;
t296 = t248 * t249;
t298 = t246 * t251;
t189 = t229 * t296 - t298;
t228 = sin(t234);
t300 = t228 * t248;
t330 = -t334 * t188 - t336 * t189 - t332 * t300;
t190 = -t229 * t298 + t296;
t191 = t229 * t295 + t297;
t299 = t228 * t251;
t329 = -t334 * t190 - t336 * t191 - t332 * t299;
t328 = t335 * t188 + t337 * t189 - t334 * t300;
t327 = t335 * t190 + t337 * t191 - t334 * t299;
t326 = t337 * t188 + t338 * t189 - t336 * t300;
t325 = t337 * t190 + t338 * t191 - t336 * t299;
t324 = t332 * t229 + (t334 * t246 - t336 * t249) * t228;
t323 = t334 * t229 + (-t335 * t246 + t337 * t249) * t228;
t322 = t336 * t229 + (-t337 * t246 + t338 * t249) * t228;
t303 = Icges(4,4) * t233;
t273 = -Icges(4,2) * t232 + t303;
t170 = -Icges(4,6) * t251 + t248 * t273;
t171 = Icges(4,6) * t248 + t251 * t273;
t304 = Icges(4,4) * t232;
t276 = Icges(4,1) * t233 - t304;
t172 = -Icges(4,5) * t251 + t248 * t276;
t173 = Icges(4,5) * t248 + t251 * t276;
t305 = Icges(3,4) * t250;
t274 = -Icges(3,2) * t247 + t305;
t182 = -Icges(3,6) * t251 + t248 * t274;
t183 = Icges(3,6) * t248 + t251 * t274;
t306 = Icges(3,4) * t247;
t277 = Icges(3,1) * t250 - t306;
t184 = -Icges(3,5) * t251 + t248 * t277;
t185 = Icges(3,5) * t248 + t251 * t277;
t200 = Icges(4,2) * t233 + t304;
t201 = Icges(4,1) * t232 + t303;
t215 = Icges(3,2) * t250 + t306;
t218 = Icges(3,1) * t247 + t305;
t226 = -qJD(2) * t251 + V_base(5);
t227 = qJD(2) * t248 + V_base(4);
t235 = V_base(6) + qJD(1);
t321 = (-t200 * t232 + t201 * t233 - t215 * t247 + t218 * t250) * t235 + (-t171 * t232 + t173 * t233 - t183 * t247 + t185 * t250) * t227 + (-t170 * t232 + t172 * t233 - t182 * t247 + t184 * t250) * t226;
t320 = (Icges(3,5) * t247 + Icges(4,5) * t232 + Icges(3,6) * t250 + Icges(4,6) * t233) * t235 + (t333 * t248 + t331 * t251) * t227 + (t331 * t248 - t333 * t251) * t226;
t313 = pkin(2) * t247;
t312 = pkin(3) * t232;
t311 = t250 * pkin(2);
t310 = pkin(5) * t249;
t307 = Icges(2,4) * t248;
t302 = Icges(5,4) * t228;
t301 = Icges(5,4) * t229;
t264 = qJ(6) * t228 + t229 * t310;
t294 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t300 - pkin(5) * t298 + t248 * t264;
t293 = t191 * rSges(7,1) + t190 * rSges(7,2) + rSges(7,3) * t299 + pkin(5) * t297 + t251 * t264;
t292 = (-qJ(6) - rSges(7,3)) * t229 + (rSges(7,1) * t249 - rSges(7,2) * t246 + t310) * t228;
t166 = -qJ(3) * t251 + t248 * t311;
t224 = t248 * pkin(1) - pkin(7) * t251;
t291 = -t166 - t224;
t290 = pkin(3) * t233;
t288 = qJD(5) * t228;
t287 = qJD(6) * t228;
t286 = V_base(5) * pkin(6) + V_base(1);
t139 = -pkin(8) * t251 + t248 * t290;
t283 = -t139 + t291;
t205 = qJD(4) * t248 + t227;
t282 = qJD(3) * t248 + t226 * t313 + t286;
t281 = pkin(4) * t229 + pkin(9) * t228;
t280 = rSges(3,1) * t250 - rSges(3,2) * t247;
t279 = rSges(4,1) * t233 - rSges(4,2) * t232;
t278 = rSges(5,1) * t229 - rSges(5,2) * t228;
t275 = Icges(5,1) * t229 - t302;
t272 = -Icges(5,2) * t228 + t301;
t269 = Icges(5,5) * t229 - Icges(5,6) * t228;
t268 = t226 * t312 + t282;
t225 = pkin(1) * t251 + t248 * pkin(7);
t267 = -V_base(4) * pkin(6) + t235 * t225 + V_base(2);
t266 = V_base(4) * t224 - t225 * V_base(5) + V_base(3);
t204 = V_base(5) + (-qJD(2) - qJD(4)) * t251;
t265 = t227 * t166 + t266;
t263 = (-Icges(5,3) * t251 + t248 * t269) * t204 + (Icges(5,3) * t248 + t251 * t269) * t205 + (Icges(5,5) * t228 + Icges(5,6) * t229) * t235;
t167 = qJ(3) * t248 + t251 * t311;
t260 = -qJD(3) * t251 + t235 * t167 + t267;
t140 = pkin(8) * t248 + t251 * t290;
t259 = t227 * t139 + (-t140 - t167) * t226 + t265;
t178 = t281 * t248;
t196 = pkin(4) * t228 - pkin(9) * t229;
t258 = t204 * t196 + (-t178 + t283) * t235 + t268;
t179 = t281 * t251;
t257 = t205 * t178 - t179 * t204 + t259;
t256 = t235 * t140 + (-t312 - t313) * t227 + t260;
t255 = t235 * t179 - t205 * t196 + t256;
t158 = -Icges(5,6) * t251 + t248 * t272;
t159 = Icges(5,6) * t248 + t251 * t272;
t160 = -Icges(5,5) * t251 + t248 * t275;
t161 = Icges(5,5) * t248 + t251 * t275;
t193 = Icges(5,2) * t229 + t302;
t194 = Icges(5,1) * t228 + t301;
t254 = (-t159 * t228 + t161 * t229) * t205 + (-t158 * t228 + t160 * t229) * t204 + (-t193 * t228 + t194 * t229) * t235;
t239 = Icges(2,4) * t251;
t223 = rSges(2,1) * t251 - t248 * rSges(2,2);
t222 = t248 * rSges(2,1) + rSges(2,2) * t251;
t221 = rSges(3,1) * t247 + rSges(3,2) * t250;
t220 = Icges(2,1) * t251 - t307;
t219 = Icges(2,1) * t248 + t239;
t217 = -Icges(2,2) * t248 + t239;
t216 = Icges(2,2) * t251 + t307;
t211 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t210 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t209 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t203 = -qJD(5) * t229 + t235;
t202 = rSges(4,1) * t232 + rSges(4,2) * t233;
t195 = rSges(5,1) * t228 + rSges(5,2) * t229;
t187 = t248 * rSges(3,3) + t251 * t280;
t186 = -rSges(3,3) * t251 + t248 * t280;
t177 = t248 * rSges(4,3) + t251 * t279;
t176 = -rSges(4,3) * t251 + t248 * t279;
t175 = t251 * t288 + t205;
t174 = t248 * t288 + t204;
t165 = t248 * rSges(5,3) + t251 * t278;
t164 = -rSges(5,3) * t251 + t248 * t278;
t163 = V_base(5) * rSges(2,3) - t222 * t235 + t286;
t162 = t223 * t235 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = t222 * V_base(4) - t223 * V_base(5) + V_base(3);
t152 = -rSges(6,3) * t229 + (rSges(6,1) * t249 - rSges(6,2) * t246) * t228;
t136 = t191 * rSges(6,1) + t190 * rSges(6,2) + rSges(6,3) * t299;
t134 = rSges(6,1) * t189 + rSges(6,2) * t188 + rSges(6,3) * t300;
t118 = t221 * t226 + (-t186 - t224) * t235 + t286;
t117 = t187 * t235 - t221 * t227 + t267;
t116 = t186 * t227 - t187 * t226 + t266;
t115 = t202 * t226 + (-t176 + t291) * t235 + t282;
t114 = t235 * t177 + (-t202 - t313) * t227 + t260;
t113 = t176 * t227 + (-t167 - t177) * t226 + t265;
t112 = t195 * t204 + (-t164 + t283) * t235 + t268;
t111 = t235 * t165 - t205 * t195 + t256;
t110 = t164 * t205 - t165 * t204 + t259;
t109 = -t134 * t203 + t152 * t174 + t258;
t108 = t203 * t136 - t175 * t152 + t255;
t107 = t134 * t175 - t136 * t174 + t257;
t106 = t174 * t292 - t203 * t294 + t251 * t287 + t258;
t105 = -t175 * t292 + t203 * t293 + t248 * t287 + t255;
t104 = -qJD(6) * t229 - t174 * t293 + t175 * t294 + t257;
t1 = t205 * (t263 * t248 + t254 * t251) / 0.2e1 + t204 * (t254 * t248 - t263 * t251) / 0.2e1 + m(1) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(2) * (t154 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + ((t188 * t323 + t189 * t322 + t300 * t324) * t203 + (t188 * t327 + t189 * t325 + t300 * t329) * t175 + (t328 * t188 + t326 * t189 + t330 * t300) * t174) * t174 / 0.2e1 + ((t190 * t323 + t191 * t322 + t299 * t324) * t203 + (t327 * t190 + t325 * t191 + t329 * t299) * t175 + (t328 * t190 + t326 * t191 + t299 * t330) * t174) * t175 / 0.2e1 + ((-t174 * t330 - t329 * t175 - t324 * t203) * t229 + ((-t246 * t323 + t249 * t322) * t203 + (-t246 * t327 + t249 * t325) * t175 + (-t246 * t328 + t249 * t326) * t174) * t228) * t203 / 0.2e1 + (t321 * t248 - t320 * t251) * t226 / 0.2e1 + (t320 * t248 + t321 * t251) * t227 / 0.2e1 + ((-t248 * t216 + t219 * t251 + Icges(1,4)) * V_base(5) + (-t248 * t217 + t220 * t251 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t216 * t251 + t248 * t219 + Icges(1,2)) * V_base(5) + (t217 * t251 + t248 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t159 * t229 + t161 * t228) * t205 + (t158 * t229 + t160 * t228) * t204 + (t171 * t233 + t173 * t232 + t183 * t250 + t185 * t247) * t227 + (t170 * t233 + t172 * t232 + t182 * t250 + t184 * t247) * t226 + (t193 * t229 + t194 * t228 + t200 * t233 + t201 * t232 + t215 * t250 + t218 * t247 + Icges(2,3)) * t235) * t235 / 0.2e1 + t235 * V_base(4) * (Icges(2,5) * t251 - Icges(2,6) * t248) + t235 * V_base(5) * (Icges(2,5) * t248 + Icges(2,6) * t251) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
