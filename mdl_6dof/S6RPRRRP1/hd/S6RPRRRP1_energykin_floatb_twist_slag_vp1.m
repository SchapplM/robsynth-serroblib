% Calculate kinetic energy for
% S6RPRRRP1
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:21
% EndTime: 2019-03-09 05:55:23
% DurationCPUTime: 2.46s
% Computational Cost: add. (2160->293), mult. (1766->412), div. (0->0), fcn. (1612->10), ass. (0->151)
t305 = Icges(6,1) + Icges(7,1);
t304 = -Icges(6,4) + Icges(7,5);
t303 = Icges(7,4) + Icges(6,5);
t302 = Icges(6,2) + Icges(7,3);
t301 = -Icges(7,6) + Icges(6,6);
t300 = -Icges(6,3) - Icges(7,2);
t299 = rSges(7,1) + pkin(5);
t298 = rSges(7,3) + qJ(6);
t225 = qJ(1) + pkin(10);
t217 = sin(t225);
t218 = cos(t225);
t230 = cos(qJ(5));
t226 = qJ(3) + qJ(4);
t221 = cos(t226);
t227 = sin(qJ(5));
t268 = t221 * t227;
t167 = t217 * t268 + t218 * t230;
t267 = t221 * t230;
t168 = t217 * t267 - t218 * t227;
t220 = sin(t226);
t270 = t217 * t220;
t297 = t302 * t167 + t304 * t168 - t301 * t270;
t169 = -t217 * t230 + t218 * t268;
t170 = t217 * t227 + t218 * t267;
t269 = t218 * t220;
t296 = t302 * t169 + t304 * t170 - t301 * t269;
t295 = -t301 * t167 + t303 * t168 - t300 * t270;
t294 = -t301 * t169 + t303 * t170 - t300 * t269;
t293 = t304 * t167 + t305 * t168 + t303 * t270;
t292 = t304 * t169 + t305 * t170 + t303 * t269;
t291 = t301 * t221 + (t302 * t227 + t304 * t230) * t220;
t290 = t300 * t221 + (-t301 * t227 + t303 * t230) * t220;
t289 = -t303 * t221 + (t304 * t227 + t305 * t230) * t220;
t229 = sin(qJ(1));
t282 = pkin(1) * t229;
t232 = cos(qJ(1));
t281 = pkin(1) * t232;
t228 = sin(qJ(3));
t280 = pkin(3) * t228;
t231 = cos(qJ(3));
t279 = pkin(3) * t231;
t278 = -pkin(6) - qJ(2);
t276 = Icges(2,4) * t229;
t275 = Icges(3,4) * t217;
t274 = Icges(4,4) * t228;
t273 = Icges(4,4) * t231;
t272 = Icges(5,4) * t220;
t271 = Icges(5,4) * t221;
t266 = rSges(7,2) * t270 + t298 * t167 + t299 * t168;
t265 = rSges(7,2) * t269 + t298 * t169 + t299 * t170;
t264 = -rSges(7,2) * t221 + (t298 * t227 + t299 * t230) * t220;
t263 = qJD(5) * t220;
t219 = V_base(6) + qJD(1);
t262 = t219 * t281 + V_base(2);
t261 = V_base(5) * pkin(6) + V_base(1);
t197 = qJD(3) * t217 + V_base(4);
t185 = t217 * pkin(2) - t218 * pkin(7);
t258 = -t185 - t282;
t257 = V_base(5) * qJ(2) + t261;
t256 = V_base(4) * t282 + qJD(2) + V_base(3);
t174 = qJD(4) * t217 + t197;
t132 = -pkin(8) * t218 + t217 * t279;
t255 = -t132 + t258;
t196 = -qJD(3) * t218 + V_base(5);
t254 = t196 * t280 + t257;
t253 = pkin(4) * t221 + pkin(9) * t220;
t252 = rSges(4,1) * t231 - rSges(4,2) * t228;
t251 = rSges(5,1) * t221 - rSges(5,2) * t220;
t250 = Icges(4,1) * t231 - t274;
t249 = Icges(5,1) * t221 - t272;
t248 = -Icges(4,2) * t228 + t273;
t247 = -Icges(5,2) * t220 + t271;
t246 = Icges(4,5) * t231 - Icges(4,6) * t228;
t245 = Icges(5,5) * t221 - Icges(5,6) * t220;
t173 = V_base(5) + (-qJD(3) - qJD(4)) * t218;
t244 = (-Icges(5,3) * t218 + t217 * t245) * t173 + (Icges(5,3) * t217 + t218 * t245) * t174 + (Icges(5,5) * t220 + Icges(5,6) * t221) * t219;
t243 = (-Icges(4,3) * t218 + t217 * t246) * t196 + (Icges(4,3) * t217 + t218 * t246) * t197 + (Icges(4,5) * t228 + Icges(4,6) * t231) * t219;
t186 = t218 * pkin(2) + t217 * pkin(7);
t242 = t219 * t186 + t278 * V_base(4) + t262;
t241 = V_base(4) * t185 + (-t186 - t281) * V_base(5) + t256;
t133 = pkin(8) * t217 + t218 * t279;
t240 = t219 * t133 - t197 * t280 + t242;
t165 = t253 * t217;
t191 = pkin(4) * t220 - pkin(9) * t221;
t239 = t173 * t191 + (-t165 + t255) * t219 + t254;
t238 = t197 * t132 - t133 * t196 + t241;
t166 = t253 * t218;
t237 = t219 * t166 - t174 * t191 + t240;
t236 = t174 * t165 - t166 * t173 + t238;
t137 = -Icges(5,6) * t218 + t217 * t247;
t138 = Icges(5,6) * t217 + t218 * t247;
t139 = -Icges(5,5) * t218 + t217 * t249;
t140 = Icges(5,5) * t217 + t218 * t249;
t188 = Icges(5,2) * t221 + t272;
t189 = Icges(5,1) * t220 + t271;
t235 = (-t138 * t220 + t140 * t221) * t174 + (-t137 * t220 + t139 * t221) * t173 + (-t188 * t220 + t189 * t221) * t219;
t149 = -Icges(4,6) * t218 + t217 * t248;
t150 = Icges(4,6) * t217 + t218 * t248;
t151 = -Icges(4,5) * t218 + t217 * t250;
t152 = Icges(4,5) * t217 + t218 * t250;
t201 = Icges(4,2) * t231 + t274;
t204 = Icges(4,1) * t228 + t273;
t234 = (-t150 * t228 + t152 * t231) * t197 + (-t149 * t228 + t151 * t231) * t196 + (-t201 * t228 + t204 * t231) * t219;
t223 = Icges(2,4) * t232;
t215 = Icges(3,4) * t218;
t209 = rSges(2,1) * t232 - rSges(2,2) * t229;
t208 = rSges(2,1) * t229 + rSges(2,2) * t232;
t207 = rSges(4,1) * t228 + rSges(4,2) * t231;
t206 = Icges(2,1) * t232 - t276;
t205 = Icges(2,1) * t229 + t223;
t203 = -Icges(2,2) * t229 + t223;
t202 = Icges(2,2) * t232 + t276;
t195 = -qJD(5) * t221 + t219;
t194 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t193 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t192 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t190 = rSges(5,1) * t220 + rSges(5,2) * t221;
t184 = rSges(3,1) * t218 - rSges(3,2) * t217;
t183 = rSges(3,1) * t217 + rSges(3,2) * t218;
t182 = Icges(3,1) * t218 - t275;
t181 = Icges(3,1) * t217 + t215;
t180 = -Icges(3,2) * t217 + t215;
t179 = Icges(3,2) * t218 + t275;
t164 = -rSges(6,3) * t221 + (rSges(6,1) * t230 - rSges(6,2) * t227) * t220;
t156 = rSges(4,3) * t217 + t218 * t252;
t155 = -rSges(4,3) * t218 + t217 * t252;
t154 = V_base(5) * rSges(2,3) - t208 * t219 + t261;
t153 = t209 * t219 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t208 * V_base(4) - t209 * V_base(5) + V_base(3);
t144 = t218 * t263 + t174;
t143 = t217 * t263 + t173;
t142 = rSges(5,3) * t217 + t218 * t251;
t141 = -rSges(5,3) * t218 + t217 * t251;
t128 = V_base(5) * rSges(3,3) + (-t183 - t282) * t219 + t257;
t127 = t184 * t219 + (-rSges(3,3) + t278) * V_base(4) + t262;
t124 = t183 * V_base(4) + (-t184 - t281) * V_base(5) + t256;
t123 = rSges(6,1) * t170 - rSges(6,2) * t169 + rSges(6,3) * t269;
t121 = rSges(6,1) * t168 - rSges(6,2) * t167 + rSges(6,3) * t270;
t107 = t196 * t207 + (-t155 + t258) * t219 + t257;
t106 = t156 * t219 - t197 * t207 + t242;
t105 = t155 * t197 - t156 * t196 + t241;
t104 = t173 * t190 + (-t141 + t255) * t219 + t254;
t103 = t142 * t219 - t174 * t190 + t240;
t102 = t141 * t174 - t142 * t173 + t238;
t101 = -t121 * t195 + t143 * t164 + t239;
t100 = t123 * t195 - t144 * t164 + t237;
t99 = t121 * t144 - t123 * t143 + t236;
t98 = qJD(6) * t169 + t143 * t264 - t195 * t266 + t239;
t97 = qJD(6) * t167 - t144 * t264 + t195 * t265 + t237;
t96 = qJD(6) * t220 * t227 - t143 * t265 + t144 * t266 + t236;
t1 = t197 * (t243 * t217 + t234 * t218) / 0.2e1 + t196 * (t234 * t217 - t243 * t218) / 0.2e1 + t174 * (t244 * t217 + t235 * t218) / 0.2e1 + t173 * (t235 * t217 - t244 * t218) / 0.2e1 + m(1) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(2) * (t145 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(3) * (t124 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + ((t167 * t291 + t168 * t289 + t270 * t290) * t195 + (t167 * t296 + t168 * t292 + t270 * t294) * t144 + (t297 * t167 + t293 * t168 + t295 * t270) * t143) * t143 / 0.2e1 + ((t169 * t291 + t170 * t289 + t269 * t290) * t195 + (t296 * t169 + t292 * t170 + t294 * t269) * t144 + (t169 * t297 + t293 * t170 + t295 * t269) * t143) * t144 / 0.2e1 + ((-t143 * t295 - t144 * t294 - t195 * t290) * t221 + ((t227 * t291 + t230 * t289) * t195 + (t227 * t296 + t230 * t292) * t144 + (t227 * t297 + t293 * t230) * t143) * t220) * t195 / 0.2e1 + ((-t179 * t217 + t181 * t218 - t202 * t229 + t205 * t232 + Icges(1,4)) * V_base(5) + (-t217 * t180 + t218 * t182 - t229 * t203 + t232 * t206 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t218 * t179 + t217 * t181 + t232 * t202 + t229 * t205 + Icges(1,2)) * V_base(5) + (t180 * t218 + t182 * t217 + t203 * t232 + t206 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t150 * t231 + t152 * t228) * t197 + (t149 * t231 + t151 * t228) * t196 + (t138 * t221 + t140 * t220) * t174 + (t137 * t221 + t139 * t220) * t173 + (t221 * t188 + t220 * t189 + t231 * t201 + t228 * t204 + Icges(2,3) + Icges(3,3)) * t219) * t219 / 0.2e1 + t219 * V_base(5) * (Icges(2,5) * t229 + Icges(3,5) * t217 + Icges(2,6) * t232 + Icges(3,6) * t218) + t219 * V_base(4) * (Icges(2,5) * t232 + Icges(3,5) * t218 - Icges(2,6) * t229 - Icges(3,6) * t217) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
