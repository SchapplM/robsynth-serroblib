% Calculate kinetic energy for
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:19
% EndTime: 2019-03-09 12:45:23
% DurationCPUTime: 4.51s
% Computational Cost: add. (1550->315), mult. (2327->454), div. (0->0), fcn. (2221->8), ass. (0->157)
t355 = Icges(3,4) + Icges(4,6);
t354 = Icges(3,1) + Icges(4,2);
t353 = -Icges(3,2) - Icges(4,3);
t258 = cos(qJ(2));
t352 = t355 * t258;
t255 = sin(qJ(2));
t351 = t355 * t255;
t350 = Icges(4,4) - Icges(3,5);
t349 = Icges(4,5) - Icges(3,6);
t348 = t255 * t353 + t352;
t347 = t258 * t354 - t351;
t346 = Icges(4,1) + Icges(3,3);
t345 = Icges(6,1) + Icges(7,1);
t344 = Icges(6,4) + Icges(7,4);
t343 = Icges(7,5) + Icges(6,5);
t342 = Icges(6,2) + Icges(7,2);
t341 = Icges(7,6) + Icges(6,6);
t340 = Icges(7,3) + Icges(6,3);
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t339 = t348 * t256 + t349 * t259;
t338 = -t349 * t256 + t348 * t259;
t337 = t347 * t256 + t350 * t259;
t336 = -t350 * t256 + t347 * t259;
t335 = t258 * t353 - t351;
t334 = t255 * t354 + t352;
t333 = t349 * t255 - t350 * t258;
t253 = qJ(4) + qJ(5);
t247 = sin(t253);
t248 = cos(t253);
t305 = t255 * t259;
t191 = -t247 * t256 + t248 * t305;
t192 = t247 * t305 + t248 * t256;
t301 = t258 * t259;
t332 = t341 * t191 + t343 * t192 + t340 * t301;
t306 = t255 * t256;
t193 = t247 * t259 + t248 * t306;
t194 = t247 * t306 - t248 * t259;
t303 = t256 * t258;
t331 = t341 * t193 + t343 * t194 + t340 * t303;
t330 = t342 * t191 + t344 * t192 + t341 * t301;
t329 = t342 * t193 + t344 * t194 + t341 * t303;
t328 = t344 * t191 + t345 * t192 + t343 * t301;
t327 = t344 * t193 + t345 * t194 + t343 * t303;
t326 = (-t343 * t247 - t341 * t248) * t258 + t340 * t255;
t325 = (-t344 * t247 - t342 * t248) * t258 + t341 * t255;
t324 = (-t345 * t247 - t344 * t248) * t258 + t343 * t255;
t237 = -qJD(2) * t259 + V_base(5);
t238 = qJD(2) * t256 + V_base(4);
t244 = V_base(6) + qJD(1);
t323 = (t335 * t255 + t334 * t258) * t244 + (-t338 * t255 + t336 * t258) * t238 + (-t339 * t255 + t337 * t258) * t237;
t322 = (-t350 * t255 - t349 * t258) * t244 + (t346 * t256 + t333 * t259) * t238 + (t333 * t256 - t346 * t259) * t237;
t254 = sin(qJ(4));
t315 = pkin(4) * t254;
t314 = t255 * pkin(8);
t257 = cos(qJ(4));
t313 = t257 * pkin(4);
t311 = Icges(2,4) * t256;
t304 = t256 * t257;
t302 = t257 * t259;
t287 = pkin(5) * t247;
t270 = qJ(6) * t258 + t255 * t287;
t296 = pkin(5) * t248;
t300 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t301 + t256 * t296 + t259 * t270;
t299 = rSges(7,1) * t194 + rSges(7,2) * t193 + rSges(7,3) * t303 + t256 * t270 - t259 * t296;
t298 = (-rSges(7,1) * t247 - rSges(7,2) * t248 - t287) * t258 + (qJ(6) + rSges(7,3)) * t255;
t283 = pkin(2) * t258 + qJ(3) * t255;
t204 = t283 * t256;
t235 = pkin(1) * t256 - pkin(7) * t259;
t297 = -t204 - t235;
t294 = qJD(3) * t255;
t293 = qJD(4) * t258;
t292 = qJD(5) * t258;
t291 = qJD(6) * t258;
t290 = V_base(5) * pkin(6) + V_base(1);
t199 = t259 * t293 + t238;
t229 = qJD(4) * t255 + t244;
t230 = pkin(2) * t255 - qJ(3) * t258;
t286 = t237 * t230 + t259 * t294 + t290;
t285 = rSges(3,1) * t258 - rSges(3,2) * t255;
t284 = -rSges(4,2) * t258 + rSges(4,3) * t255;
t198 = t256 * t293 + t237;
t236 = pkin(1) * t259 + pkin(7) * t256;
t276 = -V_base(4) * pkin(6) + t244 * t236 + V_base(2);
t275 = V_base(4) * t235 - t236 * V_base(5) + V_base(3);
t274 = pkin(9) * t258 + t255 * t315;
t205 = t283 * t259;
t271 = t244 * t205 + t256 * t294 + t276;
t269 = -qJD(3) * t258 + t238 * t204 + t275;
t212 = -t259 * pkin(3) + pkin(8) * t303;
t268 = t237 * t314 + (-t212 + t297) * t244 + t286;
t152 = t256 * t274 - t259 * t313;
t196 = pkin(9) * t255 - t258 * t315;
t267 = -t152 * t229 + t198 * t196 + t268;
t211 = t256 * pkin(3) + pkin(8) * t301;
t266 = t244 * t211 + (-t230 - t314) * t238 + t271;
t265 = t238 * t212 + (-t205 - t211) * t237 + t269;
t151 = t256 * t313 + t259 * t274;
t264 = t229 * t151 - t196 * t199 + t266;
t263 = -t151 * t198 + t199 * t152 + t265;
t249 = Icges(2,4) * t259;
t234 = rSges(2,1) * t259 - rSges(2,2) * t256;
t233 = rSges(2,1) * t256 + rSges(2,2) * t259;
t232 = rSges(3,1) * t255 + rSges(3,2) * t258;
t231 = -rSges(4,2) * t255 - rSges(4,3) * t258;
t228 = Icges(2,1) * t259 - t311;
t227 = Icges(2,1) * t256 + t249;
t225 = -Icges(2,2) * t256 + t249;
t224 = Icges(2,2) * t259 + t311;
t215 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t214 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t213 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t207 = qJD(5) * t255 + t229;
t203 = t254 * t306 - t302;
t202 = t254 * t259 + t255 * t304;
t201 = t254 * t305 + t304;
t200 = -t254 * t256 + t255 * t302;
t190 = -rSges(4,1) * t259 + t256 * t284;
t189 = rSges(4,1) * t256 + t259 * t284;
t188 = rSges(3,3) * t256 + t259 * t285;
t187 = rSges(5,3) * t255 + (-rSges(5,1) * t254 - rSges(5,2) * t257) * t258;
t186 = -rSges(3,3) * t259 + t256 * t285;
t177 = Icges(5,5) * t255 + (-Icges(5,1) * t254 - Icges(5,4) * t257) * t258;
t174 = Icges(5,6) * t255 + (-Icges(5,4) * t254 - Icges(5,2) * t257) * t258;
t171 = Icges(5,3) * t255 + (-Icges(5,5) * t254 - Icges(5,6) * t257) * t258;
t168 = t259 * t292 + t199;
t167 = t256 * t292 + t198;
t166 = rSges(6,3) * t255 + (-rSges(6,1) * t247 - rSges(6,2) * t248) * t258;
t157 = V_base(5) * rSges(2,3) - t233 * t244 + t290;
t156 = t234 * t244 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t233 * V_base(4) - t234 * V_base(5) + V_base(3);
t150 = rSges(5,1) * t203 + rSges(5,2) * t202 + rSges(5,3) * t303;
t149 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t301;
t148 = Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t303;
t147 = Icges(5,1) * t201 + Icges(5,4) * t200 + Icges(5,5) * t301;
t146 = Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t303;
t145 = Icges(5,4) * t201 + Icges(5,2) * t200 + Icges(5,6) * t301;
t144 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t303;
t143 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t301;
t141 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t303;
t139 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t301;
t124 = t232 * t237 + (-t186 - t235) * t244 + t290;
t123 = t188 * t244 - t232 * t238 + t276;
t120 = t186 * t238 - t188 * t237 + t275;
t119 = t231 * t237 + (-t190 + t297) * t244 + t286;
t118 = t189 * t244 + (-t230 - t231) * t238 + t271;
t117 = t190 * t238 + (-t189 - t205) * t237 + t269;
t116 = -t150 * t229 + t187 * t198 + t268;
t115 = t149 * t229 - t187 * t199 + t266;
t114 = -t149 * t198 + t150 * t199 + t265;
t113 = -t141 * t207 + t166 * t167 + t267;
t112 = t139 * t207 - t166 * t168 + t264;
t111 = -t139 * t167 + t141 * t168 + t263;
t110 = t167 * t298 - t207 * t299 + t259 * t291 + t267;
t109 = -t168 * t298 + t207 * t300 + t256 * t291 + t264;
t108 = qJD(6) * t255 - t167 * t300 + t168 * t299 + t263;
t1 = t199 * ((t143 * t301 + t145 * t200 + t147 * t201) * t199 + (t144 * t301 + t146 * t200 + t148 * t201) * t198 + (t171 * t301 + t174 * t200 + t177 * t201) * t229) / 0.2e1 + t198 * ((t143 * t303 + t145 * t202 + t147 * t203) * t199 + (t144 * t303 + t146 * t202 + t148 * t203) * t198 + (t171 * t303 + t174 * t202 + t177 * t203) * t229) / 0.2e1 + t229 * ((t143 * t199 + t144 * t198 + t171 * t229) * t255 + ((-t145 * t257 - t147 * t254) * t199 + (-t146 * t257 - t148 * t254) * t198 + (-t174 * t257 - t177 * t254) * t229) * t258) / 0.2e1 + m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + ((t193 * t325 + t194 * t324 + t303 * t326) * t207 + (t330 * t193 + t328 * t194 + t303 * t332) * t168 + (t329 * t193 + t327 * t194 + t331 * t303) * t167) * t167 / 0.2e1 + ((t191 * t325 + t192 * t324 + t301 * t326) * t207 + (t330 * t191 + t328 * t192 + t332 * t301) * t168 + (t191 * t329 + t192 * t327 + t301 * t331) * t167) * t168 / 0.2e1 + (((-t247 * t324 - t248 * t325) * t207 + (-t247 * t328 - t248 * t330) * t168 + (-t247 * t327 - t248 * t329) * t167) * t258 + (t331 * t167 + t168 * t332 + t326 * t207) * t255) * t207 / 0.2e1 + (t323 * t256 - t322 * t259) * t237 / 0.2e1 + (t322 * t256 + t323 * t259) * t238 / 0.2e1 + ((-t224 * t256 + t227 * t259 + Icges(1,4)) * V_base(5) + (-t225 * t256 + t228 * t259 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t224 * t259 + t227 * t256 + Icges(1,2)) * V_base(5) + (t225 * t259 + t228 * t256 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t336 * t255 + t338 * t258) * t238 + (t337 * t255 + t339 * t258) * t237 + (t334 * t255 - t335 * t258 + Icges(2,3)) * t244) * t244 / 0.2e1 + t244 * V_base(4) * (Icges(2,5) * t259 - Icges(2,6) * t256) + V_base(5) * t244 * (Icges(2,5) * t256 + Icges(2,6) * t259) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
