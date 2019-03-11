% Calculate kinetic energy for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:45
% EndTime: 2019-03-09 09:58:49
% DurationCPUTime: 3.72s
% Computational Cost: add. (1483->307), mult. (2234->412), div. (0->0), fcn. (2138->8), ass. (0->153)
t355 = Icges(3,4) + Icges(4,6);
t354 = Icges(3,1) + Icges(4,2);
t353 = -Icges(3,2) - Icges(4,3);
t258 = cos(qJ(2));
t352 = t355 * t258;
t255 = sin(qJ(2));
t351 = t355 * t255;
t350 = Icges(4,4) - Icges(3,5);
t349 = Icges(4,5) - Icges(3,6);
t348 = t353 * t255 + t352;
t347 = t354 * t258 - t351;
t346 = Icges(4,1) + Icges(3,3);
t345 = Icges(6,1) + Icges(7,1);
t344 = Icges(6,4) - Icges(7,5);
t343 = Icges(7,4) + Icges(6,5);
t342 = Icges(6,2) + Icges(7,3);
t341 = Icges(7,6) - Icges(6,6);
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t340 = t256 * t348 + t259 * t349;
t339 = -t256 * t349 + t259 * t348;
t338 = t347 * t256 + t259 * t350;
t337 = -t256 * t350 + t347 * t259;
t336 = t353 * t258 - t351;
t335 = t354 * t255 + t352;
t334 = t349 * t255 - t258 * t350;
t333 = Icges(6,3) + Icges(7,2) + Icges(5,3);
t332 = rSges(7,1) + pkin(5);
t331 = rSges(7,3) + qJ(6);
t254 = sin(qJ(4));
t312 = pkin(4) * t254;
t272 = qJ(5) * t258 + t255 * t312;
t257 = cos(qJ(4));
t310 = pkin(4) * t257;
t157 = t256 * t272 - t259 * t310;
t239 = qJD(2) * t256 + V_base(4);
t289 = qJD(4) * t258;
t202 = t259 * t289 + t239;
t330 = qJD(5) * t255 + t202 * t157;
t252 = qJ(4) + pkin(9);
t246 = cos(t252);
t302 = t255 * t259;
t245 = sin(t252);
t303 = t245 * t256;
t173 = -t246 * t302 + t303;
t301 = t256 * t246;
t174 = t245 * t302 + t301;
t296 = t258 * t259;
t329 = t173 * t342 - t174 * t344 + t296 * t341;
t175 = t245 * t259 + t255 * t301;
t176 = -t246 * t259 + t255 * t303;
t298 = t256 * t258;
t328 = -t175 * t342 - t176 * t344 + t298 * t341;
t327 = -t173 * t344 + t174 * t345 + t296 * t343;
t326 = t175 * t344 + t176 * t345 + t298 * t343;
t325 = (t245 * t344 + t246 * t342) * t258 + t341 * t255;
t324 = (-t245 * t345 - t246 * t344) * t258 + t343 * t255;
t238 = -qJD(2) * t259 + V_base(5);
t247 = V_base(6) + qJD(1);
t323 = (t255 * t336 + t258 * t335) * t247 + (-t255 * t339 + t258 * t337) * t239 + (-t255 * t340 + t258 * t338) * t238;
t322 = (-t255 * t350 - t349 * t258) * t247 + (t256 * t346 + t334 * t259) * t239 + (t334 * t256 - t259 * t346) * t238;
t297 = t257 * t259;
t300 = t256 * t254;
t203 = t255 * t297 - t300;
t299 = t256 * t257;
t204 = t254 * t302 + t299;
t321 = Icges(5,5) * t204 + Icges(5,6) * t203 + t173 * t341 + t174 * t343 + t296 * t333;
t205 = t254 * t259 + t255 * t299;
t206 = t255 * t300 - t297;
t320 = Icges(5,5) * t206 + Icges(5,6) * t205 - t175 * t341 + t176 * t343 + t298 * t333;
t319 = (-Icges(5,5) * t254 - Icges(5,6) * t257 - t245 * t343 + t246 * t341) * t258 + t333 * t255;
t311 = pkin(8) * t255;
t308 = Icges(2,4) * t256;
t295 = rSges(7,2) * t296 + t331 * t173 + t332 * t174;
t294 = rSges(7,2) * t298 - t331 * t175 + t332 * t176;
t293 = rSges(7,2) * t255 + (-t332 * t245 + t331 * t246) * t258;
t281 = pkin(2) * t258 + qJ(3) * t255;
t207 = t281 * t256;
t236 = t256 * pkin(1) - pkin(7) * t259;
t292 = -t207 - t236;
t291 = qJD(3) * t255;
t290 = qJD(3) * t258;
t288 = qJD(5) * t258;
t287 = V_base(5) * pkin(6) + V_base(1);
t231 = pkin(2) * t255 - qJ(3) * t258;
t284 = t238 * t231 + t259 * t291 + t287;
t283 = rSges(3,1) * t258 - rSges(3,2) * t255;
t282 = -rSges(4,2) * t258 + rSges(4,3) * t255;
t237 = pkin(1) * t259 + t256 * pkin(7);
t274 = -V_base(4) * pkin(6) + t247 * t237 + V_base(2);
t273 = V_base(4) * t236 - t237 * V_base(5) + V_base(3);
t271 = t239 * t207 + t273;
t208 = t281 * t259;
t268 = t247 * t208 + t256 * t291 + t274;
t213 = -pkin(3) * t259 + pkin(8) * t298;
t267 = t238 * t311 + (-t213 + t292) * t247 + t284;
t212 = t256 * pkin(3) + pkin(8) * t296;
t266 = t239 * t213 + (-t208 - t212) * t238 + t271;
t199 = qJ(5) * t255 - t258 * t312;
t201 = t256 * t289 + t238;
t265 = t201 * t199 + t259 * t288 + t267;
t264 = t247 * t212 + (-t231 - t311) * t239 + t268;
t263 = t266 - t290;
t156 = t256 * t310 + t259 * t272;
t230 = qJD(4) * t255 + t247;
t262 = t230 * t156 + t256 * t288 + t264;
t250 = Icges(2,4) * t259;
t235 = rSges(2,1) * t259 - t256 * rSges(2,2);
t234 = t256 * rSges(2,1) + rSges(2,2) * t259;
t233 = rSges(3,1) * t255 + rSges(3,2) * t258;
t232 = -rSges(4,2) * t255 - rSges(4,3) * t258;
t229 = Icges(2,1) * t259 - t308;
t228 = Icges(2,1) * t256 + t250;
t226 = -Icges(2,2) * t256 + t250;
t225 = Icges(2,2) * t259 + t308;
t216 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t215 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t214 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t197 = -rSges(4,1) * t259 + t256 * t282;
t196 = t256 * rSges(4,1) + t259 * t282;
t195 = t256 * rSges(3,3) + t259 * t283;
t194 = rSges(5,3) * t255 + (-rSges(5,1) * t254 - rSges(5,2) * t257) * t258;
t193 = -rSges(3,3) * t259 + t256 * t283;
t183 = Icges(5,5) * t255 + (-Icges(5,1) * t254 - Icges(5,4) * t257) * t258;
t180 = Icges(5,6) * t255 + (-Icges(5,4) * t254 - Icges(5,2) * t257) * t258;
t169 = rSges(6,3) * t255 + (-rSges(6,1) * t245 - rSges(6,2) * t246) * t258;
t161 = V_base(5) * rSges(2,3) - t234 * t247 + t287;
t160 = t235 * t247 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t159 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t155 = rSges(5,1) * t206 + rSges(5,2) * t205 + rSges(5,3) * t298;
t154 = t204 * rSges(5,1) + t203 * rSges(5,2) + rSges(5,3) * t296;
t153 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t298;
t152 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t296;
t151 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t298;
t150 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t296;
t144 = rSges(6,1) * t176 + rSges(6,2) * t175 + rSges(6,3) * t298;
t142 = t174 * rSges(6,1) - t173 * rSges(6,2) + rSges(6,3) * t296;
t127 = t233 * t238 + (-t193 - t236) * t247 + t287;
t126 = t195 * t247 - t233 * t239 + t274;
t125 = t193 * t239 - t195 * t238 + t273;
t124 = t232 * t238 + (-t197 + t292) * t247 + t284;
t123 = t196 * t247 + (-t231 - t232) * t239 + t268;
t122 = -t290 + t197 * t239 + (-t196 - t208) * t238 + t271;
t121 = -t155 * t230 + t194 * t201 + t267;
t120 = t154 * t230 - t194 * t202 + t264;
t119 = -t154 * t201 + t155 * t202 + t263;
t118 = t169 * t201 + (-t144 - t157) * t230 + t265;
t117 = t142 * t230 + (-t169 - t199) * t202 + t262;
t116 = t144 * t202 + (-t142 - t156) * t201 + t263 + t330;
t115 = qJD(6) * t173 + t293 * t201 + (-t157 - t294) * t230 + t265;
t114 = -qJD(6) * t175 + t295 * t230 + (-t199 - t293) * t202 + t262;
t113 = (qJD(6) * t246 - qJD(3)) * t258 + t294 * t202 + (-t156 - t295) * t201 + t266 + t330;
t1 = m(2) * (t159 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(3) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(7) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(6) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(1) * (t214 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 + (t323 * t256 - t322 * t259) * t238 / 0.2e1 + (t322 * t256 + t323 * t259) * t239 / 0.2e1 + ((-t256 * t225 + t228 * t259 + Icges(1,4)) * V_base(5) + (-t256 * t226 + t259 * t229 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t259 * t225 + t256 * t228 + Icges(1,2)) * V_base(5) + (t226 * t259 + t256 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t175 * t325 + t176 * t324 + t180 * t205 + t183 * t206 + t298 * t319) * t230 + (t150 * t205 + t152 * t206 - t175 * t329 + t176 * t327 + t298 * t321) * t202 + (t205 * t151 + t206 * t153 - t328 * t175 + t326 * t176 + t320 * t298) * t201) * t201 / 0.2e1 + ((t173 * t325 + t174 * t324 + t203 * t180 + t204 * t183 + t296 * t319) * t230 + (t203 * t150 + t204 * t152 + t329 * t173 + t327 * t174 + t321 * t296) * t202 + (t203 * t151 + t204 * t153 + t173 * t328 + t174 * t326 + t296 * t320) * t201) * t202 / 0.2e1 + (((-t180 * t257 - t183 * t254 - t245 * t324 + t246 * t325) * t230 + (-t150 * t257 - t152 * t254 - t245 * t327 + t246 * t329) * t202 + (-t151 * t257 - t153 * t254 - t245 * t326 + t246 * t328) * t201) * t258 + (t201 * t320 + t202 * t321 + t230 * t319) * t255) * t230 / 0.2e1 + ((t255 * t337 + t258 * t339) * t239 + (t255 * t338 + t258 * t340) * t238 + (t335 * t255 - t336 * t258 + Icges(2,3)) * t247) * t247 / 0.2e1 + t247 * V_base(4) * (Icges(2,5) * t259 - Icges(2,6) * t256) + t247 * V_base(5) * (Icges(2,5) * t256 + Icges(2,6) * t259) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
