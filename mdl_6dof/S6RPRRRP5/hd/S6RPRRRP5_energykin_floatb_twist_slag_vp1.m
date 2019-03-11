% Calculate kinetic energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:41
% EndTime: 2019-03-09 06:10:44
% DurationCPUTime: 2.89s
% Computational Cost: add. (2107->312), mult. (1990->446), div. (0->0), fcn. (1838->10), ass. (0->166)
t332 = Icges(6,1) + Icges(7,1);
t331 = -Icges(6,4) + Icges(7,5);
t330 = Icges(7,4) + Icges(6,5);
t329 = Icges(6,2) + Icges(7,3);
t328 = -Icges(7,6) + Icges(6,6);
t327 = -Icges(6,3) - Icges(7,2);
t326 = rSges(7,1) + pkin(5);
t325 = rSges(7,3) + qJ(6);
t245 = pkin(10) + qJ(3);
t236 = qJ(4) + t245;
t231 = cos(t236);
t251 = cos(qJ(5));
t252 = cos(qJ(1));
t296 = t251 * t252;
t249 = sin(qJ(5));
t250 = sin(qJ(1));
t298 = t250 * t249;
t190 = t231 * t298 + t296;
t297 = t250 * t251;
t299 = t249 * t252;
t191 = t231 * t297 - t299;
t230 = sin(t236);
t301 = t230 * t250;
t324 = t329 * t190 + t331 * t191 - t328 * t301;
t192 = t231 * t299 - t297;
t193 = t231 * t296 + t298;
t300 = t230 * t252;
t323 = t329 * t192 + t331 * t193 - t328 * t300;
t322 = -t328 * t190 + t330 * t191 - t327 * t301;
t321 = -t328 * t192 + t330 * t193 - t327 * t300;
t320 = t331 * t190 + t332 * t191 + t330 * t301;
t319 = t331 * t192 + t332 * t193 + t330 * t300;
t318 = t328 * t231 + (t329 * t249 + t331 * t251) * t230;
t317 = t327 * t231 + (-t328 * t249 + t330 * t251) * t230;
t316 = -t330 * t231 + (t331 * t249 + t332 * t251) * t230;
t246 = sin(pkin(10));
t311 = pkin(2) * t246;
t234 = sin(t245);
t310 = pkin(3) * t234;
t247 = cos(pkin(10));
t309 = t247 * pkin(2);
t308 = Icges(2,4) * t250;
t307 = Icges(3,4) * t246;
t306 = Icges(3,4) * t247;
t305 = Icges(4,4) * t234;
t235 = cos(t245);
t304 = Icges(4,4) * t235;
t303 = Icges(5,4) * t230;
t302 = Icges(5,4) * t231;
t294 = rSges(7,2) * t301 + t325 * t190 + t191 * t326;
t293 = rSges(7,2) * t300 + t325 * t192 + t193 * t326;
t292 = -rSges(7,2) * t231 + (t325 * t249 + t251 * t326) * t230;
t167 = -pkin(7) * t252 + t250 * t309;
t223 = t250 * pkin(1) - qJ(2) * t252;
t291 = -t167 - t223;
t290 = pkin(3) * t235;
t288 = qJD(5) * t230;
t287 = V_base(4) * t223 + V_base(3);
t286 = V_base(5) * pkin(6) + V_base(1);
t141 = -pkin(8) * t252 + t250 * t290;
t283 = -t141 + t291;
t228 = qJD(3) * t250 + V_base(4);
t237 = V_base(6) + qJD(1);
t282 = qJD(2) * t250 + t286;
t207 = qJD(4) * t250 + t228;
t281 = V_base(5) * t311 + t282;
t280 = pkin(4) * t231 + pkin(9) * t230;
t279 = rSges(3,1) * t247 - rSges(3,2) * t246;
t278 = rSges(4,1) * t235 - rSges(4,2) * t234;
t277 = rSges(5,1) * t231 - rSges(5,2) * t230;
t276 = Icges(3,1) * t247 - t307;
t275 = Icges(4,1) * t235 - t305;
t274 = Icges(5,1) * t231 - t303;
t273 = -Icges(3,2) * t246 + t306;
t272 = -Icges(4,2) * t234 + t304;
t271 = -Icges(5,2) * t230 + t302;
t270 = Icges(3,5) * t247 - Icges(3,6) * t246;
t269 = Icges(4,5) * t235 - Icges(4,6) * t234;
t268 = Icges(5,5) * t231 - Icges(5,6) * t230;
t225 = pkin(1) * t252 + t250 * qJ(2);
t267 = -qJD(2) * t252 + t237 * t225 + V_base(2);
t227 = -qJD(3) * t252 + V_base(5);
t266 = t227 * t310 + t281;
t206 = V_base(5) + (-qJD(3) - qJD(4)) * t252;
t265 = (-Icges(5,3) * t252 + t250 * t268) * t206 + (Icges(5,3) * t250 + t252 * t268) * t207 + (Icges(5,5) * t230 + Icges(5,6) * t231) * t237;
t264 = (-Icges(4,3) * t252 + t250 * t269) * t227 + (Icges(4,3) * t250 + t252 * t269) * t228 + (Icges(4,5) * t234 + Icges(4,6) * t235) * t237;
t168 = pkin(7) * t250 + t252 * t309;
t263 = V_base(4) * t167 + (-t168 - t225) * V_base(5) + t287;
t262 = (-Icges(3,3) * t252 + t250 * t270) * V_base(5) + (Icges(3,3) * t250 + t252 * t270) * V_base(4) + (Icges(3,5) * t246 + Icges(3,6) * t247) * t237;
t142 = pkin(8) * t250 + t252 * t290;
t261 = t228 * t141 - t142 * t227 + t263;
t260 = t237 * t168 + (-pkin(6) - t311) * V_base(4) + t267;
t180 = t280 * t250;
t198 = pkin(4) * t230 - pkin(9) * t231;
t259 = t206 * t198 + (-t180 + t283) * t237 + t266;
t181 = t280 * t252;
t258 = t207 * t180 - t181 * t206 + t261;
t257 = t237 * t142 - t228 * t310 + t260;
t256 = t237 * t181 - t207 * t198 + t257;
t159 = -Icges(5,6) * t252 + t250 * t271;
t160 = Icges(5,6) * t250 + t252 * t271;
t161 = -Icges(5,5) * t252 + t250 * t274;
t162 = Icges(5,5) * t250 + t252 * t274;
t195 = Icges(5,2) * t231 + t303;
t196 = Icges(5,1) * t230 + t302;
t255 = (-t160 * t230 + t162 * t231) * t207 + (-t159 * t230 + t161 * t231) * t206 + (-t195 * t230 + t196 * t231) * t237;
t171 = -Icges(4,6) * t252 + t250 * t272;
t172 = Icges(4,6) * t250 + t252 * t272;
t173 = -Icges(4,5) * t252 + t250 * t275;
t174 = Icges(4,5) * t250 + t252 * t275;
t202 = Icges(4,2) * t235 + t305;
t203 = Icges(4,1) * t234 + t304;
t254 = (-t172 * t234 + t174 * t235) * t228 + (-t171 * t234 + t173 * t235) * t227 + (-t202 * t234 + t203 * t235) * t237;
t184 = -Icges(3,6) * t252 + t250 * t273;
t185 = Icges(3,6) * t250 + t252 * t273;
t186 = -Icges(3,5) * t252 + t250 * t276;
t187 = Icges(3,5) * t250 + t252 * t276;
t214 = Icges(3,2) * t247 + t307;
t215 = Icges(3,1) * t246 + t306;
t253 = (-t185 * t246 + t187 * t247) * V_base(4) + (-t184 * t246 + t186 * t247) * V_base(5) + (-t214 * t246 + t215 * t247) * t237;
t242 = Icges(2,4) * t252;
t226 = rSges(2,1) * t252 - t250 * rSges(2,2);
t224 = t250 * rSges(2,1) + rSges(2,2) * t252;
t222 = Icges(2,1) * t252 - t308;
t221 = Icges(2,1) * t250 + t242;
t220 = -Icges(2,2) * t250 + t242;
t219 = Icges(2,2) * t252 + t308;
t218 = Icges(2,5) * t252 - Icges(2,6) * t250;
t217 = Icges(2,5) * t250 + Icges(2,6) * t252;
t216 = rSges(3,1) * t246 + rSges(3,2) * t247;
t212 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t211 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t210 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t205 = -qJD(5) * t231 + t237;
t204 = rSges(4,1) * t234 + rSges(4,2) * t235;
t197 = rSges(5,1) * t230 + rSges(5,2) * t231;
t189 = t250 * rSges(3,3) + t252 * t279;
t188 = -rSges(3,3) * t252 + t250 * t279;
t178 = t250 * rSges(4,3) + t252 * t278;
t177 = -rSges(4,3) * t252 + t250 * t278;
t176 = t252 * t288 + t207;
t175 = t250 * t288 + t206;
t166 = t250 * rSges(5,3) + t252 * t277;
t165 = -rSges(5,3) * t252 + t250 * t277;
t164 = V_base(5) * rSges(2,3) - t224 * t237 + t286;
t163 = t226 * t237 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t224 * V_base(4) - t226 * V_base(5) + V_base(3);
t152 = -rSges(6,3) * t231 + (rSges(6,1) * t251 - rSges(6,2) * t249) * t230;
t136 = t193 * rSges(6,1) - t192 * rSges(6,2) + rSges(6,3) * t300;
t134 = rSges(6,1) * t191 - rSges(6,2) * t190 + rSges(6,3) * t301;
t120 = t216 * V_base(5) + (-t188 - t223) * t237 + t282;
t119 = t237 * t189 + (-pkin(6) - t216) * V_base(4) + t267;
t118 = t188 * V_base(4) + (-t189 - t225) * V_base(5) + t287;
t117 = t204 * t227 + (-t177 + t291) * t237 + t281;
t116 = t237 * t178 - t228 * t204 + t260;
t115 = t177 * t228 - t178 * t227 + t263;
t114 = t197 * t206 + (-t165 + t283) * t237 + t266;
t113 = t237 * t166 - t207 * t197 + t257;
t112 = t165 * t207 - t166 * t206 + t261;
t111 = -t134 * t205 + t152 * t175 + t259;
t110 = t205 * t136 - t176 * t152 + t256;
t109 = t134 * t176 - t136 * t175 + t258;
t108 = qJD(6) * t192 + t175 * t292 - t205 * t294 + t259;
t107 = qJD(6) * t190 - t176 * t292 + t205 * t293 + t256;
t106 = qJD(6) * t230 * t249 - t175 * t293 + t176 * t294 + t258;
t1 = m(1) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(3) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + t228 * (t264 * t250 + t254 * t252) / 0.2e1 + t227 * (t254 * t250 - t264 * t252) / 0.2e1 + t207 * (t250 * t265 + t252 * t255) / 0.2e1 + t206 * (t255 * t250 - t265 * t252) / 0.2e1 + ((t190 * t318 + t191 * t316 + t301 * t317) * t205 + (t190 * t323 + t191 * t319 + t301 * t321) * t176 + (t324 * t190 + t320 * t191 + t322 * t301) * t175) * t175 / 0.2e1 + ((t192 * t318 + t193 * t316 + t300 * t317) * t205 + (t323 * t192 + t319 * t193 + t321 * t300) * t176 + (t192 * t324 + t320 * t193 + t322 * t300) * t175) * t176 / 0.2e1 + ((-t175 * t322 - t176 * t321 - t205 * t317) * t231 + ((t249 * t318 + t251 * t316) * t205 + (t249 * t323 + t251 * t319) * t176 + (t249 * t324 + t320 * t251) * t175) * t230) * t205 / 0.2e1 + (t218 * t237 + t250 * t262 + t252 * t253 + (-t250 * t219 + t221 * t252 + Icges(1,4)) * V_base(5) + (-t250 * t220 + t222 * t252 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t217 * t237 + t250 * t253 - t252 * t262 + (t219 * t252 + t250 * t221 + Icges(1,2)) * V_base(5) + (t220 * t252 + t250 * t222 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t172 * t235 + t174 * t234) * t228 + (t171 * t235 + t173 * t234) * t227 + (t160 * t231 + t162 * t230) * t207 + (t159 * t231 + t161 * t230) * t206 + (t184 * t247 + t186 * t246 + t217) * V_base(5) + (t185 * t247 + t187 * t246 + t218) * V_base(4) + (t195 * t231 + t196 * t230 + t202 * t235 + t203 * t234 + t214 * t247 + t215 * t246 + Icges(2,3)) * t237) * t237 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
