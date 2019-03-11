% Calculate kinetic energy for
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:53
% EndTime: 2019-03-09 06:06:56
% DurationCPUTime: 3.01s
% Computational Cost: add. (2149->317), mult. (1999->453), div. (0->0), fcn. (1837->10), ass. (0->167)
t333 = Icges(6,1) + Icges(7,1);
t332 = Icges(6,4) + Icges(7,4);
t331 = -Icges(7,5) - Icges(6,5);
t330 = Icges(6,2) + Icges(7,2);
t329 = -Icges(7,6) - Icges(6,6);
t328 = -Icges(7,3) - Icges(6,3);
t243 = pkin(10) + qJ(3);
t234 = qJ(4) + t243;
t229 = cos(t234);
t250 = cos(qJ(5));
t251 = cos(qJ(1));
t297 = t250 * t251;
t248 = sin(qJ(5));
t249 = sin(qJ(1));
t299 = t249 * t248;
t188 = -t229 * t299 - t297;
t298 = t249 * t250;
t300 = t248 * t251;
t189 = t229 * t298 - t300;
t228 = sin(t234);
t302 = t228 * t249;
t327 = -t329 * t188 - t331 * t189 - t328 * t302;
t190 = -t229 * t300 + t298;
t191 = t229 * t297 + t299;
t301 = t228 * t251;
t326 = -t329 * t190 - t331 * t191 - t328 * t301;
t325 = t330 * t188 + t332 * t189 - t329 * t302;
t324 = t330 * t190 + t332 * t191 - t329 * t301;
t323 = t332 * t188 + t333 * t189 - t331 * t302;
t322 = t332 * t190 + t333 * t191 - t331 * t301;
t321 = t328 * t229 + (t329 * t248 - t331 * t250) * t228;
t320 = t329 * t229 + (-t330 * t248 + t332 * t250) * t228;
t319 = t331 * t229 + (-t332 * t248 + t333 * t250) * t228;
t244 = sin(pkin(10));
t314 = pkin(2) * t244;
t232 = sin(t243);
t313 = pkin(3) * t232;
t245 = cos(pkin(10));
t312 = t245 * pkin(2);
t311 = pkin(5) * t250;
t309 = Icges(2,4) * t249;
t308 = Icges(3,4) * t244;
t307 = Icges(3,4) * t245;
t306 = Icges(4,4) * t232;
t233 = cos(t243);
t305 = Icges(4,4) * t233;
t304 = Icges(5,4) * t228;
t303 = Icges(5,4) * t229;
t265 = qJ(6) * t228 + t229 * t311;
t295 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t302 - pkin(5) * t300 + t249 * t265;
t294 = t191 * rSges(7,1) + t190 * rSges(7,2) + rSges(7,3) * t301 + pkin(5) * t299 + t251 * t265;
t293 = (-qJ(6) - rSges(7,3)) * t229 + (rSges(7,1) * t250 - rSges(7,2) * t248 + t311) * t228;
t166 = -pkin(7) * t251 + t249 * t312;
t221 = t249 * pkin(1) - qJ(2) * t251;
t292 = -t166 - t221;
t291 = pkin(3) * t233;
t289 = qJD(5) * t228;
t288 = qJD(6) * t228;
t287 = V_base(4) * t221 + V_base(3);
t286 = V_base(5) * pkin(6) + V_base(1);
t139 = -pkin(8) * t251 + t249 * t291;
t283 = -t139 + t292;
t226 = qJD(3) * t249 + V_base(4);
t235 = V_base(6) + qJD(1);
t282 = qJD(2) * t249 + t286;
t205 = qJD(4) * t249 + t226;
t281 = V_base(5) * t314 + t282;
t280 = pkin(4) * t229 + pkin(9) * t228;
t279 = rSges(3,1) * t245 - rSges(3,2) * t244;
t278 = rSges(4,1) * t233 - rSges(4,2) * t232;
t277 = rSges(5,1) * t229 - rSges(5,2) * t228;
t276 = Icges(3,1) * t245 - t308;
t275 = Icges(4,1) * t233 - t306;
t274 = Icges(5,1) * t229 - t304;
t273 = -Icges(3,2) * t244 + t307;
t272 = -Icges(4,2) * t232 + t305;
t271 = -Icges(5,2) * t228 + t303;
t270 = Icges(3,5) * t245 - Icges(3,6) * t244;
t269 = Icges(4,5) * t233 - Icges(4,6) * t232;
t268 = Icges(5,5) * t229 - Icges(5,6) * t228;
t223 = pkin(1) * t251 + t249 * qJ(2);
t267 = -qJD(2) * t251 + t235 * t223 + V_base(2);
t225 = -qJD(3) * t251 + V_base(5);
t266 = t225 * t313 + t281;
t204 = V_base(5) + (-qJD(3) - qJD(4)) * t251;
t264 = (-Icges(5,3) * t251 + t249 * t268) * t204 + (Icges(5,3) * t249 + t251 * t268) * t205 + (Icges(5,5) * t228 + Icges(5,6) * t229) * t235;
t263 = (-Icges(4,3) * t251 + t249 * t269) * t225 + (Icges(4,3) * t249 + t251 * t269) * t226 + (Icges(4,5) * t232 + Icges(4,6) * t233) * t235;
t167 = pkin(7) * t249 + t251 * t312;
t262 = V_base(4) * t166 + (-t167 - t223) * V_base(5) + t287;
t261 = (-Icges(3,3) * t251 + t249 * t270) * V_base(5) + (Icges(3,3) * t249 + t251 * t270) * V_base(4) + (Icges(3,5) * t244 + Icges(3,6) * t245) * t235;
t140 = pkin(8) * t249 + t251 * t291;
t260 = t226 * t139 - t140 * t225 + t262;
t259 = t235 * t167 + (-pkin(6) - t314) * V_base(4) + t267;
t178 = t280 * t249;
t196 = pkin(4) * t228 - pkin(9) * t229;
t258 = t204 * t196 + (-t178 + t283) * t235 + t266;
t179 = t280 * t251;
t257 = t205 * t178 - t179 * t204 + t260;
t256 = t235 * t140 - t226 * t313 + t259;
t255 = t235 * t179 - t205 * t196 + t256;
t158 = -Icges(5,6) * t251 + t249 * t271;
t159 = Icges(5,6) * t249 + t251 * t271;
t160 = -Icges(5,5) * t251 + t249 * t274;
t161 = Icges(5,5) * t249 + t251 * t274;
t193 = Icges(5,2) * t229 + t304;
t194 = Icges(5,1) * t228 + t303;
t254 = (-t159 * t228 + t161 * t229) * t205 + (-t158 * t228 + t160 * t229) * t204 + (-t193 * t228 + t194 * t229) * t235;
t170 = -Icges(4,6) * t251 + t249 * t272;
t171 = Icges(4,6) * t249 + t251 * t272;
t172 = -Icges(4,5) * t251 + t249 * t275;
t173 = Icges(4,5) * t249 + t251 * t275;
t200 = Icges(4,2) * t233 + t306;
t201 = Icges(4,1) * t232 + t305;
t253 = (-t171 * t232 + t173 * t233) * t226 + (-t170 * t232 + t172 * t233) * t225 + (-t200 * t232 + t201 * t233) * t235;
t182 = -Icges(3,6) * t251 + t249 * t273;
t183 = Icges(3,6) * t249 + t251 * t273;
t184 = -Icges(3,5) * t251 + t249 * t276;
t185 = Icges(3,5) * t249 + t251 * t276;
t212 = Icges(3,2) * t245 + t308;
t213 = Icges(3,1) * t244 + t307;
t252 = (-t183 * t244 + t185 * t245) * V_base(4) + (-t182 * t244 + t184 * t245) * V_base(5) + (-t212 * t244 + t213 * t245) * t235;
t240 = Icges(2,4) * t251;
t224 = rSges(2,1) * t251 - t249 * rSges(2,2);
t222 = t249 * rSges(2,1) + rSges(2,2) * t251;
t220 = Icges(2,1) * t251 - t309;
t219 = Icges(2,1) * t249 + t240;
t218 = -Icges(2,2) * t249 + t240;
t217 = Icges(2,2) * t251 + t309;
t216 = Icges(2,5) * t251 - Icges(2,6) * t249;
t215 = Icges(2,5) * t249 + Icges(2,6) * t251;
t214 = rSges(3,1) * t244 + rSges(3,2) * t245;
t210 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t209 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t208 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t203 = -qJD(5) * t229 + t235;
t202 = rSges(4,1) * t232 + rSges(4,2) * t233;
t195 = rSges(5,1) * t228 + rSges(5,2) * t229;
t187 = t249 * rSges(3,3) + t251 * t279;
t186 = -rSges(3,3) * t251 + t249 * t279;
t177 = t249 * rSges(4,3) + t251 * t278;
t176 = -rSges(4,3) * t251 + t249 * t278;
t175 = t251 * t289 + t205;
t174 = t249 * t289 + t204;
t165 = t249 * rSges(5,3) + t251 * t277;
t164 = -rSges(5,3) * t251 + t249 * t277;
t163 = V_base(5) * rSges(2,3) - t222 * t235 + t286;
t162 = t224 * t235 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = t222 * V_base(4) - t224 * V_base(5) + V_base(3);
t151 = -rSges(6,3) * t229 + (rSges(6,1) * t250 - rSges(6,2) * t248) * t228;
t136 = t191 * rSges(6,1) + t190 * rSges(6,2) + rSges(6,3) * t301;
t134 = rSges(6,1) * t189 + rSges(6,2) * t188 + rSges(6,3) * t302;
t118 = t214 * V_base(5) + (-t186 - t221) * t235 + t282;
t117 = t235 * t187 + (-pkin(6) - t214) * V_base(4) + t267;
t116 = t186 * V_base(4) + (-t187 - t223) * V_base(5) + t287;
t115 = t202 * t225 + (-t176 + t292) * t235 + t281;
t114 = t235 * t177 - t226 * t202 + t259;
t113 = t176 * t226 - t177 * t225 + t262;
t112 = t195 * t204 + (-t164 + t283) * t235 + t266;
t111 = t235 * t165 - t205 * t195 + t256;
t110 = t164 * t205 - t165 * t204 + t260;
t109 = -t134 * t203 + t151 * t174 + t258;
t108 = t203 * t136 - t175 * t151 + t255;
t107 = t134 * t175 - t136 * t174 + t257;
t106 = t174 * t293 - t203 * t295 + t251 * t288 + t258;
t105 = -t175 * t293 + t203 * t294 + t249 * t288 + t255;
t104 = -qJD(6) * t229 - t174 * t294 + t175 * t295 + t257;
t1 = t226 * (t263 * t249 + t253 * t251) / 0.2e1 + t225 * (t253 * t249 - t263 * t251) / 0.2e1 + t205 * (t264 * t249 + t254 * t251) / 0.2e1 + t204 * (t254 * t249 - t264 * t251) / 0.2e1 + m(1) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(2) * (t154 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + ((t188 * t320 + t189 * t319 + t302 * t321) * t203 + (t188 * t324 + t189 * t322 + t302 * t326) * t175 + (t325 * t188 + t323 * t189 + t302 * t327) * t174) * t174 / 0.2e1 + ((t190 * t320 + t191 * t319 + t301 * t321) * t203 + (t190 * t324 + t191 * t322 + t301 * t326) * t175 + (t325 * t190 + t323 * t191 + t301 * t327) * t174) * t175 / 0.2e1 + ((-t327 * t174 - t326 * t175 - t321 * t203) * t229 + ((-t248 * t320 + t250 * t319) * t203 + (-t248 * t324 + t250 * t322) * t175 + (-t248 * t325 + t250 * t323) * t174) * t228) * t203 / 0.2e1 + (t216 * t235 + t261 * t249 + t252 * t251 + (-t249 * t217 + t219 * t251 + Icges(1,4)) * V_base(5) + (-t249 * t218 + t220 * t251 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t215 * t235 + t252 * t249 - t261 * t251 + (t217 * t251 + t249 * t219 + Icges(1,2)) * V_base(5) + (t218 * t251 + t249 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t171 * t233 + t173 * t232) * t226 + (t170 * t233 + t172 * t232) * t225 + (t159 * t229 + t161 * t228) * t205 + (t158 * t229 + t160 * t228) * t204 + (t182 * t245 + t184 * t244 + t215) * V_base(5) + (t183 * t245 + t185 * t244 + t216) * V_base(4) + (t193 * t229 + t194 * t228 + t200 * t233 + t201 * t232 + t212 * t245 + t213 * t244 + Icges(2,3)) * t235) * t235 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
