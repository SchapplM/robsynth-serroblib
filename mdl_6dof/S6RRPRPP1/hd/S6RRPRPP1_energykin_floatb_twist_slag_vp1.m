% Calculate kinetic energy for
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:03
% EndTime: 2019-03-09 09:45:06
% DurationCPUTime: 2.90s
% Computational Cost: add. (2182->323), mult. (2260->440), div. (0->0), fcn. (2164->10), ass. (0->159)
t340 = Icges(6,1) + Icges(7,1);
t339 = -Icges(6,4) + Icges(7,5);
t338 = Icges(7,4) + Icges(6,5);
t337 = Icges(6,2) + Icges(7,3);
t336 = -Icges(7,6) + Icges(6,6);
t335 = Icges(3,3) + Icges(4,3);
t252 = qJ(2) + pkin(9);
t243 = sin(t252);
t245 = cos(t252);
t256 = sin(qJ(2));
t259 = cos(qJ(2));
t334 = Icges(3,5) * t259 + Icges(4,5) * t245 - Icges(3,6) * t256 - Icges(4,6) * t243;
t333 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t332 = rSges(7,1) + pkin(5);
t331 = rSges(7,3) + qJ(6);
t251 = qJ(4) + pkin(10);
t244 = cos(t251);
t260 = cos(qJ(1));
t242 = sin(t251);
t257 = sin(qJ(1));
t299 = t257 * t242;
t189 = t244 * t260 + t245 * t299;
t298 = t257 * t244;
t190 = -t242 * t260 + t245 * t298;
t303 = t243 * t257;
t330 = t337 * t189 + t339 * t190 - t336 * t303;
t301 = t245 * t260;
t191 = t242 * t301 - t298;
t192 = t244 * t301 + t299;
t302 = t243 * t260;
t329 = t337 * t191 + t339 * t192 - t336 * t302;
t328 = t339 * t189 + t340 * t190 + t338 * t303;
t327 = t339 * t191 + t340 * t192 + t338 * t302;
t326 = t336 * t245 + (t337 * t242 + t339 * t244) * t243;
t325 = -t338 * t245 + (t339 * t242 + t340 * t244) * t243;
t304 = Icges(4,4) * t245;
t278 = -Icges(4,2) * t243 + t304;
t181 = -Icges(4,6) * t260 + t257 * t278;
t182 = Icges(4,6) * t257 + t260 * t278;
t305 = Icges(4,4) * t243;
t280 = Icges(4,1) * t245 - t305;
t183 = -Icges(4,5) * t260 + t257 * t280;
t184 = Icges(4,5) * t257 + t260 * t280;
t306 = Icges(3,4) * t259;
t279 = -Icges(3,2) * t256 + t306;
t195 = -Icges(3,6) * t260 + t257 * t279;
t196 = Icges(3,6) * t257 + t260 * t279;
t307 = Icges(3,4) * t256;
t281 = Icges(3,1) * t259 - t307;
t197 = -Icges(3,5) * t260 + t257 * t281;
t198 = Icges(3,5) * t257 + t260 * t281;
t211 = Icges(4,2) * t245 + t305;
t212 = Icges(4,1) * t243 + t304;
t224 = Icges(3,2) * t259 + t307;
t227 = Icges(3,1) * t256 + t306;
t238 = -qJD(2) * t260 + V_base(5);
t239 = qJD(2) * t257 + V_base(4);
t246 = V_base(6) + qJD(1);
t324 = (-t211 * t243 + t212 * t245 - t224 * t256 + t227 * t259) * t246 + (-t182 * t243 + t184 * t245 - t196 * t256 + t198 * t259) * t239 + (-t181 * t243 + t183 * t245 - t195 * t256 + t197 * t259) * t238;
t323 = (Icges(3,5) * t256 + Icges(4,5) * t243 + Icges(3,6) * t259 + Icges(4,6) * t245) * t246 + (t335 * t257 + t334 * t260) * t239 + (t334 * t257 - t335 * t260) * t238;
t258 = cos(qJ(4));
t295 = t258 * t260;
t255 = sin(qJ(4));
t297 = t257 * t255;
t205 = -t245 * t297 - t295;
t296 = t257 * t258;
t300 = t255 * t260;
t206 = t245 * t296 - t300;
t322 = Icges(5,5) * t206 + Icges(5,6) * t205 - t336 * t189 + t338 * t190 - t333 * t303;
t207 = -t245 * t300 + t296;
t208 = t245 * t295 + t297;
t321 = Icges(5,5) * t208 + Icges(5,6) * t207 - t336 * t191 + t338 * t192 - t333 * t302;
t320 = t333 * t245 + (Icges(5,5) * t258 - Icges(5,6) * t255 - t336 * t242 + t338 * t244) * t243;
t313 = pkin(2) * t256;
t312 = pkin(2) * t259;
t311 = pkin(4) * t258;
t308 = Icges(2,4) * t257;
t294 = rSges(7,2) * t303 + t331 * t189 + t332 * t190;
t293 = rSges(7,2) * t302 + t331 * t191 + t332 * t192;
t292 = -rSges(7,2) * t245 + (t331 * t242 + t332 * t244) * t243;
t176 = -qJ(3) * t260 + t257 * t312;
t235 = t257 * pkin(1) - pkin(7) * t260;
t291 = -t176 - t235;
t290 = qJD(4) * t243;
t289 = qJD(5) * t243;
t288 = V_base(5) * pkin(6) + V_base(1);
t285 = qJD(3) * t257 + t238 * t313 + t288;
t284 = pkin(3) * t245 + pkin(8) * t243;
t283 = rSges(3,1) * t259 - rSges(3,2) * t256;
t282 = rSges(4,1) * t245 - rSges(4,2) * t243;
t236 = pkin(1) * t260 + t257 * pkin(7);
t275 = -V_base(4) * pkin(6) + t246 * t236 + V_base(2);
t274 = V_base(4) * t235 - t236 * V_base(5) + V_base(3);
t273 = t239 * t176 + t274;
t272 = qJ(5) * t243 + t245 * t311;
t177 = qJ(3) * t257 + t260 * t312;
t269 = -qJD(3) * t260 + t246 * t177 + t275;
t201 = t284 * t257;
t214 = pkin(3) * t243 - pkin(8) * t245;
t268 = t238 * t214 + (-t201 + t291) * t246 + t285;
t202 = t284 * t260;
t267 = t239 * t201 + (-t177 - t202) * t238 + t273;
t157 = -qJ(5) * t245 + t243 * t311;
t203 = t257 * t290 + t238;
t266 = t203 * t157 + t260 * t289 + t268;
t265 = t246 * t202 + (-t214 - t313) * t239 + t269;
t145 = -pkin(4) * t300 + t257 * t272;
t204 = t260 * t290 + t239;
t264 = -qJD(5) * t245 + t204 * t145 + t267;
t146 = pkin(4) * t297 + t260 * t272;
t217 = -qJD(4) * t245 + t246;
t263 = t217 * t146 + t257 * t289 + t265;
t249 = Icges(2,4) * t260;
t234 = rSges(2,1) * t260 - t257 * rSges(2,2);
t233 = t257 * rSges(2,1) + rSges(2,2) * t260;
t232 = rSges(3,1) * t256 + rSges(3,2) * t259;
t229 = Icges(2,1) * t260 - t308;
t228 = Icges(2,1) * t257 + t249;
t226 = -Icges(2,2) * t257 + t249;
t225 = Icges(2,2) * t260 + t308;
t220 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t219 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t218 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = rSges(4,1) * t243 + rSges(4,2) * t245;
t200 = t257 * rSges(3,3) + t260 * t283;
t199 = -rSges(3,3) * t260 + t257 * t283;
t187 = t257 * rSges(4,3) + t260 * t282;
t186 = -rSges(4,3) * t260 + t257 * t282;
t175 = -rSges(5,3) * t245 + (rSges(5,1) * t258 - rSges(5,2) * t255) * t243;
t174 = V_base(5) * rSges(2,3) - t233 * t246 + t288;
t173 = t234 * t246 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t172 = -Icges(5,5) * t245 + (Icges(5,1) * t258 - Icges(5,4) * t255) * t243;
t171 = -Icges(5,6) * t245 + (Icges(5,4) * t258 - Icges(5,2) * t255) * t243;
t168 = t233 * V_base(4) - t234 * V_base(5) + V_base(3);
t166 = -rSges(6,3) * t245 + (rSges(6,1) * t244 - rSges(6,2) * t242) * t243;
t154 = t208 * rSges(5,1) + t207 * rSges(5,2) + rSges(5,3) * t302;
t153 = rSges(5,1) * t206 + rSges(5,2) * t205 + rSges(5,3) * t303;
t152 = Icges(5,1) * t208 + Icges(5,4) * t207 + Icges(5,5) * t302;
t151 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t303;
t150 = Icges(5,4) * t208 + Icges(5,2) * t207 + Icges(5,6) * t302;
t149 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t303;
t143 = t192 * rSges(6,1) - t191 * rSges(6,2) + rSges(6,3) * t302;
t141 = rSges(6,1) * t190 - rSges(6,2) * t189 + rSges(6,3) * t303;
t126 = t232 * t238 + (-t199 - t235) * t246 + t288;
t125 = t200 * t246 - t232 * t239 + t275;
t123 = t199 * t239 - t200 * t238 + t274;
t122 = t213 * t238 + (-t186 + t291) * t246 + t285;
t121 = t246 * t187 + (-t213 - t313) * t239 + t269;
t120 = t186 * t239 + (-t177 - t187) * t238 + t273;
t119 = -t153 * t217 + t175 * t203 + t268;
t118 = t217 * t154 - t204 * t175 + t265;
t117 = t153 * t204 - t154 * t203 + t267;
t116 = t166 * t203 + (-t141 - t145) * t217 + t266;
t115 = t217 * t143 + (-t157 - t166) * t204 + t263;
t114 = t141 * t204 + (-t143 - t146) * t203 + t264;
t113 = qJD(6) * t191 + t292 * t203 + (-t145 - t294) * t217 + t266;
t112 = qJD(6) * t189 + t293 * t217 + (-t157 - t292) * t204 + t263;
t111 = qJD(6) * t242 * t243 + t294 * t204 + (-t146 - t293) * t203 + t264;
t1 = m(1) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(6) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(4) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t168 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + (t324 * t257 - t323 * t260) * t238 / 0.2e1 + (t323 * t257 + t324 * t260) * t239 / 0.2e1 + ((-t257 * t225 + t228 * t260 + Icges(1,4)) * V_base(5) + (-t257 * t226 + t260 * t229 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t260 * t225 + t257 * t228 + Icges(1,2)) * V_base(5) + (t226 * t260 + t257 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t171 * t205 + t172 * t206 + t189 * t326 + t190 * t325 + t303 * t320) * t217 + (t150 * t205 + t152 * t206 + t189 * t329 + t190 * t327 + t303 * t321) * t204 + (t205 * t149 + t206 * t151 + t330 * t189 + t328 * t190 + t322 * t303) * t203) * t203 / 0.2e1 + ((t207 * t171 + t208 * t172 + t191 * t326 + t192 * t325 + t302 * t320) * t217 + (t207 * t150 + t208 * t152 + t329 * t191 + t327 * t192 + t321 * t302) * t204 + (t207 * t149 + t208 * t151 + t191 * t330 + t328 * t192 + t322 * t302) * t203) * t204 / 0.2e1 + ((-t203 * t322 - t204 * t321 - t320 * t217) * t245 + ((-t171 * t255 + t172 * t258 + t242 * t326 + t325 * t244) * t217 + (-t150 * t255 + t152 * t258 + t242 * t329 + t244 * t327) * t204 + (-t149 * t255 + t151 * t258 + t242 * t330 + t328 * t244) * t203) * t243) * t217 / 0.2e1 + ((t182 * t245 + t184 * t243 + t196 * t259 + t198 * t256) * t239 + (t181 * t245 + t183 * t243 + t195 * t259 + t197 * t256) * t238 + (t245 * t211 + t243 * t212 + t259 * t224 + t256 * t227 + Icges(2,3)) * t246) * t246 / 0.2e1 + t246 * V_base(4) * (Icges(2,5) * t260 - Icges(2,6) * t257) + V_base(5) * t246 * (Icges(2,5) * t257 + Icges(2,6) * t260) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
