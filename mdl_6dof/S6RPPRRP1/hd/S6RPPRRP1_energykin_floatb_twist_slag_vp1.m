% Calculate kinetic energy for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:33
% EndTime: 2019-03-09 01:57:35
% DurationCPUTime: 2.48s
% Computational Cost: add. (2064->301), mult. (1711->409), div. (0->0), fcn. (1547->10), ass. (0->156)
t312 = Icges(6,1) + Icges(7,1);
t311 = Icges(6,4) + Icges(7,4);
t310 = -Icges(7,5) - Icges(6,5);
t309 = Icges(6,2) + Icges(7,2);
t308 = -Icges(7,6) - Icges(6,6);
t307 = -Icges(7,3) - Icges(6,3);
t225 = pkin(10) + qJ(4);
t219 = cos(t225);
t226 = qJ(1) + pkin(9);
t220 = cos(t226);
t233 = cos(qJ(5));
t273 = t220 * t233;
t218 = sin(t226);
t231 = sin(qJ(5));
t276 = t218 * t231;
t170 = -t219 * t276 - t273;
t274 = t220 * t231;
t275 = t218 * t233;
t171 = t219 * t275 - t274;
t217 = sin(t225);
t278 = t217 * t218;
t306 = -t308 * t170 - t310 * t171 - t307 * t278;
t172 = -t219 * t274 + t275;
t173 = t219 * t273 + t276;
t277 = t217 * t220;
t305 = -t308 * t172 - t310 * t173 - t307 * t277;
t304 = t309 * t170 + t311 * t171 - t308 * t278;
t303 = t309 * t172 + t311 * t173 - t308 * t277;
t302 = t311 * t170 + t312 * t171 - t310 * t278;
t301 = t311 * t172 + t312 * t173 - t310 * t277;
t300 = t307 * t219 + (t308 * t231 - t310 * t233) * t217;
t299 = t308 * t219 + (-t309 * t231 + t311 * t233) * t217;
t298 = t310 * t219 + (-t311 * t231 + t312 * t233) * t217;
t232 = sin(qJ(1));
t234 = cos(qJ(1));
t297 = Icges(2,5) * t232 + Icges(3,5) * t218 + Icges(2,6) * t234 + Icges(3,6) * t220;
t296 = Icges(2,5) * t234 + Icges(3,5) * t220 - Icges(2,6) * t232 - Icges(3,6) * t218;
t291 = pkin(1) * t232;
t290 = pkin(1) * t234;
t227 = sin(pkin(10));
t289 = pkin(3) * t227;
t228 = cos(pkin(10));
t288 = pkin(3) * t228;
t287 = pkin(5) * t233;
t286 = -pkin(6) - qJ(2);
t284 = Icges(2,4) * t232;
t283 = Icges(3,4) * t218;
t282 = Icges(4,4) * t227;
t281 = Icges(4,4) * t228;
t280 = Icges(5,4) * t217;
t279 = Icges(5,4) * t219;
t244 = qJ(6) * t217 + t219 * t287;
t271 = rSges(7,1) * t171 + rSges(7,2) * t170 + rSges(7,3) * t278 - pkin(5) * t274 + t218 * t244;
t270 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t277 + pkin(5) * t276 + t220 * t244;
t269 = (-qJ(6) - rSges(7,3)) * t219 + (rSges(7,1) * t233 - rSges(7,2) * t231 + t287) * t217;
t268 = qJD(5) * t217;
t267 = qJD(6) * t217;
t221 = V_base(6) + qJD(1);
t266 = t221 * t290 + V_base(2);
t265 = V_base(5) * pkin(6) + V_base(1);
t199 = qJD(4) * t218 + V_base(4);
t186 = pkin(2) * t218 - qJ(3) * t220;
t262 = -t186 - t291;
t188 = pkin(2) * t220 + qJ(3) * t218;
t261 = -t188 - t290;
t260 = V_base(5) * qJ(2) + t265;
t259 = V_base(4) * t291 + qJD(2) + V_base(3);
t134 = -pkin(7) * t220 + t218 * t288;
t258 = -t134 + t262;
t257 = qJD(3) * t218 + t260;
t256 = pkin(4) * t219 + pkin(8) * t217;
t255 = V_base(4) * t186 + t259;
t198 = -qJD(4) * t220 + V_base(5);
t254 = rSges(4,1) * t228 - rSges(4,2) * t227;
t253 = rSges(5,1) * t219 - rSges(5,2) * t217;
t252 = Icges(4,1) * t228 - t282;
t251 = Icges(5,1) * t219 - t280;
t250 = -Icges(4,2) * t227 + t281;
t249 = -Icges(5,2) * t217 + t279;
t248 = Icges(4,5) * t228 - Icges(4,6) * t227;
t247 = Icges(5,5) * t219 - Icges(5,6) * t217;
t246 = V_base(5) * t289 + t257;
t245 = -qJD(3) * t220 + t221 * t188 + t266;
t243 = (-Icges(5,3) * t220 + t218 * t247) * t198 + (Icges(5,3) * t218 + t220 * t247) * t199 + (Icges(5,5) * t217 + Icges(5,6) * t219) * t221;
t242 = (-Icges(4,3) * t220 + t218 * t248) * V_base(5) + (Icges(4,3) * t218 + t220 * t248) * V_base(4) + (Icges(4,5) * t227 + Icges(4,6) * t228) * t221;
t135 = pkin(7) * t218 + t220 * t288;
t241 = V_base(4) * t134 + (-t135 + t261) * V_base(5) + t255;
t166 = t256 * t218;
t190 = pkin(4) * t217 - pkin(8) * t219;
t240 = t198 * t190 + (-t166 + t258) * t221 + t246;
t239 = t221 * t135 + (t286 - t289) * V_base(4) + t245;
t167 = t256 * t220;
t238 = t199 * t166 - t198 * t167 + t241;
t237 = t221 * t167 - t190 * t199 + t239;
t139 = -Icges(5,6) * t220 + t218 * t249;
t140 = Icges(5,6) * t218 + t220 * t249;
t141 = -Icges(5,5) * t220 + t218 * t251;
t142 = Icges(5,5) * t218 + t220 * t251;
t179 = Icges(5,2) * t219 + t280;
t182 = Icges(5,1) * t217 + t279;
t236 = (-t140 * t217 + t142 * t219) * t199 + (-t139 * t217 + t141 * t219) * t198 + (-t179 * t217 + t182 * t219) * t221;
t149 = -Icges(4,6) * t220 + t218 * t250;
t150 = Icges(4,6) * t218 + t220 * t250;
t151 = -Icges(4,5) * t220 + t218 * t252;
t152 = Icges(4,5) * t218 + t220 * t252;
t196 = Icges(4,2) * t228 + t282;
t197 = Icges(4,1) * t227 + t281;
t235 = (-t150 * t227 + t152 * t228) * V_base(4) + (-t149 * t227 + t151 * t228) * V_base(5) + (-t196 * t227 + t197 * t228) * t221;
t223 = Icges(2,4) * t234;
t214 = Icges(3,4) * t220;
t208 = rSges(2,1) * t234 - t232 * rSges(2,2);
t207 = t232 * rSges(2,1) + rSges(2,2) * t234;
t206 = Icges(2,1) * t234 - t284;
t205 = Icges(2,1) * t232 + t223;
t204 = -Icges(2,2) * t232 + t223;
t203 = Icges(2,2) * t234 + t284;
t200 = rSges(4,1) * t227 + rSges(4,2) * t228;
t194 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t193 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t192 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t191 = -qJD(5) * t219 + t221;
t189 = rSges(3,1) * t220 - rSges(3,2) * t218;
t187 = rSges(3,1) * t218 + rSges(3,2) * t220;
t185 = rSges(5,1) * t217 + rSges(5,2) * t219;
t184 = Icges(3,1) * t220 - t283;
t183 = Icges(3,1) * t218 + t214;
t181 = -Icges(3,2) * t218 + t214;
t180 = Icges(3,2) * t220 + t283;
t169 = t220 * t268 + t199;
t168 = t218 * t268 + t198;
t164 = -rSges(6,3) * t219 + (rSges(6,1) * t233 - rSges(6,2) * t231) * t217;
t162 = V_base(5) * rSges(2,3) - t207 * t221 + t265;
t161 = t208 * t221 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = rSges(4,3) * t218 + t220 * t254;
t153 = -rSges(4,3) * t220 + t218 * t254;
t146 = t207 * V_base(4) - t208 * V_base(5) + V_base(3);
t144 = rSges(5,3) * t218 + t220 * t253;
t143 = -rSges(5,3) * t220 + t218 * t253;
t130 = V_base(5) * rSges(3,3) + (-t187 - t291) * t221 + t260;
t129 = t189 * t221 + (-rSges(3,3) + t286) * V_base(4) + t266;
t128 = V_base(4) * t187 + (-t189 - t290) * V_base(5) + t259;
t127 = rSges(6,1) * t173 + rSges(6,2) * t172 + rSges(6,3) * t277;
t125 = rSges(6,1) * t171 + rSges(6,2) * t170 + rSges(6,3) * t278;
t109 = t200 * V_base(5) + (-t153 + t262) * t221 + t257;
t108 = t154 * t221 + (-t200 + t286) * V_base(4) + t245;
t107 = V_base(4) * t153 + (-t154 + t261) * V_base(5) + t255;
t106 = t185 * t198 + (-t143 + t258) * t221 + t246;
t105 = t144 * t221 - t185 * t199 + t239;
t104 = t199 * t143 - t198 * t144 + t241;
t103 = -t125 * t191 + t164 * t168 + t240;
t102 = t127 * t191 - t164 * t169 + t237;
t101 = t169 * t125 - t168 * t127 + t238;
t100 = t168 * t269 - t191 * t271 + t220 * t267 + t240;
t99 = -t169 * t269 + t191 * t270 + t218 * t267 + t237;
t98 = -qJD(6) * t219 - t168 * t270 + t169 * t271 + t238;
t1 = t199 * (t243 * t218 + t236 * t220) / 0.2e1 + t198 * (t236 * t218 - t243 * t220) / 0.2e1 + m(1) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(2) * (t146 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(3) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + ((t299 * t170 + t298 * t171 + t300 * t278) * t191 + (t303 * t170 + t301 * t171 + t305 * t278) * t169 + (t304 * t170 + t302 * t171 + t306 * t278) * t168) * t168 / 0.2e1 + ((t299 * t172 + t298 * t173 + t300 * t277) * t191 + (t303 * t172 + t301 * t173 + t305 * t277) * t169 + (t304 * t172 + t302 * t173 + t306 * t277) * t168) * t169 / 0.2e1 + ((-t306 * t168 - t305 * t169 - t300 * t191) * t219 + ((-t299 * t231 + t298 * t233) * t191 + (-t303 * t231 + t301 * t233) * t169 + (-t304 * t231 + t302 * t233) * t168) * t217) * t191 / 0.2e1 + ((t140 * t219 + t142 * t217) * t199 + (t139 * t219 + t141 * t217) * t198 + (t149 * t228 + t151 * t227 + t297) * V_base(5) + (t150 * t228 + t152 * t227 + t296) * V_base(4) + (t219 * t179 + t217 * t182 + t228 * t196 + t197 * t227 + Icges(2,3) + Icges(3,3)) * t221) * t221 / 0.2e1 + (t242 * t218 + t235 * t220 + t296 * t221 + (-t180 * t218 + t183 * t220 - t232 * t203 + t205 * t234 + Icges(1,4)) * V_base(5) + (-t218 * t181 + t220 * t184 - t232 * t204 + t234 * t206 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t235 * t218 - t242 * t220 + t297 * t221 + (t220 * t180 + t218 * t183 + t234 * t203 + t232 * t205 + Icges(1,2)) * V_base(5) + (t181 * t220 + t184 * t218 + t204 * t234 + t232 * t206 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
