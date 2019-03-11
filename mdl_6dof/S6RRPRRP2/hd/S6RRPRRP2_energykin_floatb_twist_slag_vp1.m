% Calculate kinetic energy for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:00
% EndTime: 2019-03-09 11:43:03
% DurationCPUTime: 2.80s
% Computational Cost: add. (2133->306), mult. (2016->432), div. (0->0), fcn. (1864->10), ass. (0->161)
t337 = Icges(6,1) + Icges(7,1);
t336 = -Icges(6,4) + Icges(7,5);
t335 = Icges(7,4) + Icges(6,5);
t334 = Icges(6,2) + Icges(7,3);
t333 = -Icges(7,6) + Icges(6,6);
t332 = Icges(3,3) + Icges(4,3);
t331 = -Icges(6,3) - Icges(7,2);
t245 = qJ(2) + pkin(10);
t234 = sin(t245);
t235 = cos(t245);
t248 = sin(qJ(2));
t251 = cos(qJ(2));
t330 = Icges(3,5) * t251 + Icges(4,5) * t235 - Icges(3,6) * t248 - Icges(4,6) * t234;
t329 = rSges(7,1) + pkin(5);
t328 = rSges(7,3) + qJ(6);
t236 = qJ(4) + t245;
t231 = cos(t236);
t250 = cos(qJ(5));
t252 = cos(qJ(1));
t294 = t250 * t252;
t247 = sin(qJ(5));
t249 = sin(qJ(1));
t296 = t249 * t247;
t190 = t231 * t296 + t294;
t295 = t249 * t250;
t297 = t247 * t252;
t191 = t231 * t295 - t297;
t230 = sin(t236);
t299 = t230 * t249;
t327 = t334 * t190 + t336 * t191 - t333 * t299;
t192 = t231 * t297 - t295;
t193 = t231 * t294 + t296;
t298 = t230 * t252;
t326 = t334 * t192 + t336 * t193 - t333 * t298;
t325 = -t333 * t190 + t335 * t191 - t331 * t299;
t324 = -t333 * t192 + t335 * t193 - t331 * t298;
t323 = t336 * t190 + t337 * t191 + t335 * t299;
t322 = t336 * t192 + t337 * t193 + t335 * t298;
t321 = t333 * t231 + (t334 * t247 + t336 * t250) * t230;
t320 = t331 * t231 + (-t333 * t247 + t335 * t250) * t230;
t319 = -t335 * t231 + (t336 * t247 + t337 * t250) * t230;
t302 = Icges(4,4) * t235;
t273 = -Icges(4,2) * t234 + t302;
t171 = -Icges(4,6) * t252 + t249 * t273;
t172 = Icges(4,6) * t249 + t252 * t273;
t303 = Icges(4,4) * t234;
t276 = Icges(4,1) * t235 - t303;
t173 = -Icges(4,5) * t252 + t249 * t276;
t174 = Icges(4,5) * t249 + t252 * t276;
t304 = Icges(3,4) * t251;
t274 = -Icges(3,2) * t248 + t304;
t184 = -Icges(3,6) * t252 + t249 * t274;
t185 = Icges(3,6) * t249 + t252 * t274;
t305 = Icges(3,4) * t248;
t277 = Icges(3,1) * t251 - t305;
t186 = -Icges(3,5) * t252 + t249 * t277;
t187 = Icges(3,5) * t249 + t252 * t277;
t202 = Icges(4,2) * t235 + t303;
t203 = Icges(4,1) * t234 + t302;
t217 = Icges(3,2) * t251 + t305;
t220 = Icges(3,1) * t248 + t304;
t228 = -qJD(2) * t252 + V_base(5);
t229 = qJD(2) * t249 + V_base(4);
t237 = V_base(6) + qJD(1);
t318 = (-t202 * t234 + t203 * t235 - t217 * t248 + t220 * t251) * t237 + (-t172 * t234 + t174 * t235 - t185 * t248 + t187 * t251) * t229 + (-t171 * t234 + t173 * t235 - t184 * t248 + t186 * t251) * t228;
t317 = (Icges(3,5) * t248 + Icges(4,5) * t234 + Icges(3,6) * t251 + Icges(4,6) * t235) * t237 + (t332 * t249 + t330 * t252) * t229 + (t330 * t249 - t332 * t252) * t228;
t310 = pkin(2) * t248;
t309 = pkin(3) * t234;
t308 = t251 * pkin(2);
t306 = Icges(2,4) * t249;
t301 = Icges(5,4) * t230;
t300 = Icges(5,4) * t231;
t293 = rSges(7,2) * t299 + t328 * t190 + t191 * t329;
t292 = rSges(7,2) * t298 + t328 * t192 + t193 * t329;
t291 = -rSges(7,2) * t231 + (t328 * t247 + t250 * t329) * t230;
t167 = -qJ(3) * t252 + t249 * t308;
t226 = t249 * pkin(1) - pkin(7) * t252;
t290 = -t167 - t226;
t289 = pkin(3) * t235;
t287 = qJD(5) * t230;
t286 = V_base(5) * pkin(6) + V_base(1);
t141 = -pkin(8) * t252 + t249 * t289;
t283 = -t141 + t290;
t207 = qJD(4) * t249 + t229;
t282 = qJD(3) * t249 + t228 * t310 + t286;
t281 = pkin(4) * t231 + pkin(9) * t230;
t280 = rSges(3,1) * t251 - rSges(3,2) * t248;
t279 = rSges(4,1) * t235 - rSges(4,2) * t234;
t278 = rSges(5,1) * t231 - rSges(5,2) * t230;
t275 = Icges(5,1) * t231 - t301;
t272 = -Icges(5,2) * t230 + t300;
t269 = Icges(5,5) * t231 - Icges(5,6) * t230;
t268 = t228 * t309 + t282;
t227 = pkin(1) * t252 + t249 * pkin(7);
t267 = -V_base(4) * pkin(6) + t237 * t227 + V_base(2);
t266 = V_base(4) * t226 - t227 * V_base(5) + V_base(3);
t206 = V_base(5) + (-qJD(2) - qJD(4)) * t252;
t265 = t229 * t167 + t266;
t264 = (-Icges(5,3) * t252 + t249 * t269) * t206 + (Icges(5,3) * t249 + t252 * t269) * t207 + (Icges(5,5) * t230 + Icges(5,6) * t231) * t237;
t168 = qJ(3) * t249 + t252 * t308;
t261 = -qJD(3) * t252 + t237 * t168 + t267;
t142 = pkin(8) * t249 + t252 * t289;
t260 = t229 * t141 + (-t142 - t168) * t228 + t265;
t180 = t281 * t249;
t198 = pkin(4) * t230 - pkin(9) * t231;
t259 = t206 * t198 + (-t180 + t283) * t237 + t268;
t181 = t281 * t252;
t258 = t207 * t180 - t181 * t206 + t260;
t257 = t237 * t142 + (-t309 - t310) * t229 + t261;
t256 = t237 * t181 - t207 * t198 + t257;
t159 = -Icges(5,6) * t252 + t249 * t272;
t160 = Icges(5,6) * t249 + t252 * t272;
t161 = -Icges(5,5) * t252 + t249 * t275;
t162 = Icges(5,5) * t249 + t252 * t275;
t195 = Icges(5,2) * t231 + t301;
t196 = Icges(5,1) * t230 + t300;
t255 = (-t160 * t230 + t162 * t231) * t207 + (-t159 * t230 + t161 * t231) * t206 + (-t195 * t230 + t196 * t231) * t237;
t241 = Icges(2,4) * t252;
t225 = rSges(2,1) * t252 - t249 * rSges(2,2);
t224 = t249 * rSges(2,1) + rSges(2,2) * t252;
t223 = rSges(3,1) * t248 + rSges(3,2) * t251;
t222 = Icges(2,1) * t252 - t306;
t221 = Icges(2,1) * t249 + t241;
t219 = -Icges(2,2) * t249 + t241;
t218 = Icges(2,2) * t252 + t306;
t213 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t212 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t211 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t205 = -qJD(5) * t231 + t237;
t204 = rSges(4,1) * t234 + rSges(4,2) * t235;
t197 = rSges(5,1) * t230 + rSges(5,2) * t231;
t189 = t249 * rSges(3,3) + t252 * t280;
t188 = -rSges(3,3) * t252 + t249 * t280;
t178 = t249 * rSges(4,3) + t252 * t279;
t177 = -rSges(4,3) * t252 + t249 * t279;
t176 = t252 * t287 + t207;
t175 = t249 * t287 + t206;
t166 = t249 * rSges(5,3) + t252 * t278;
t165 = -rSges(5,3) * t252 + t249 * t278;
t164 = V_base(5) * rSges(2,3) - t224 * t237 + t286;
t163 = t225 * t237 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t224 * V_base(4) - t225 * V_base(5) + V_base(3);
t153 = -rSges(6,3) * t231 + (rSges(6,1) * t250 - rSges(6,2) * t247) * t230;
t136 = t193 * rSges(6,1) - t192 * rSges(6,2) + rSges(6,3) * t298;
t134 = rSges(6,1) * t191 - rSges(6,2) * t190 + rSges(6,3) * t299;
t120 = t223 * t228 + (-t188 - t226) * t237 + t286;
t119 = t189 * t237 - t223 * t229 + t267;
t118 = t188 * t229 - t189 * t228 + t266;
t117 = t204 * t228 + (-t177 + t290) * t237 + t282;
t116 = t237 * t178 + (-t204 - t310) * t229 + t261;
t115 = t177 * t229 + (-t168 - t178) * t228 + t265;
t114 = t197 * t206 + (-t165 + t283) * t237 + t268;
t113 = t237 * t166 - t207 * t197 + t257;
t112 = t165 * t207 - t166 * t206 + t260;
t111 = -t134 * t205 + t153 * t175 + t259;
t110 = t205 * t136 - t176 * t153 + t256;
t109 = t134 * t176 - t136 * t175 + t258;
t108 = qJD(6) * t192 + t175 * t291 - t205 * t293 + t259;
t107 = qJD(6) * t190 - t176 * t291 + t205 * t292 + t256;
t106 = qJD(6) * t230 * t247 - t175 * t292 + t176 * t293 + t258;
t1 = m(1) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(3) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t207 * (t249 * t264 + t252 * t255) / 0.2e1 + t206 * (t255 * t249 - t264 * t252) / 0.2e1 + ((t190 * t321 + t191 * t319 + t299 * t320) * t205 + (t190 * t326 + t191 * t322 + t299 * t324) * t176 + (t327 * t190 + t323 * t191 + t325 * t299) * t175) * t175 / 0.2e1 + ((t192 * t321 + t193 * t319 + t298 * t320) * t205 + (t326 * t192 + t322 * t193 + t324 * t298) * t176 + (t192 * t327 + t323 * t193 + t325 * t298) * t175) * t176 / 0.2e1 + ((-t175 * t325 - t176 * t324 - t205 * t320) * t231 + ((t247 * t321 + t250 * t319) * t205 + (t247 * t326 + t250 * t322) * t176 + (t247 * t327 + t323 * t250) * t175) * t230) * t205 / 0.2e1 + (t318 * t249 - t317 * t252) * t228 / 0.2e1 + (t317 * t249 + t318 * t252) * t229 / 0.2e1 + ((-t249 * t218 + t221 * t252 + Icges(1,4)) * V_base(5) + (-t249 * t219 + t222 * t252 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t218 * t252 + t249 * t221 + Icges(1,2)) * V_base(5) + (t219 * t252 + t249 * t222 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t160 * t231 + t162 * t230) * t207 + (t159 * t231 + t161 * t230) * t206 + (t172 * t235 + t174 * t234 + t185 * t251 + t187 * t248) * t229 + (t171 * t235 + t173 * t234 + t184 * t251 + t186 * t248) * t228 + (t195 * t231 + t196 * t230 + t202 * t235 + t203 * t234 + t217 * t251 + t220 * t248 + Icges(2,3)) * t237) * t237 / 0.2e1 + t237 * V_base(4) * (Icges(2,5) * t252 - Icges(2,6) * t249) + V_base(5) * t237 * (Icges(2,5) * t249 + Icges(2,6) * t252) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
