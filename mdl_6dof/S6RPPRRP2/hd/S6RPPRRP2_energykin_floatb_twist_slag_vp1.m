% Calculate kinetic energy for
% S6RPPRRP2
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:01
% EndTime: 2019-03-09 02:00:04
% DurationCPUTime: 2.40s
% Computational Cost: add. (2032->296), mult. (1702->402), div. (0->0), fcn. (1548->10), ass. (0->155)
t311 = Icges(6,1) + Icges(7,1);
t310 = -Icges(6,4) + Icges(7,5);
t309 = Icges(7,4) + Icges(6,5);
t308 = Icges(6,2) + Icges(7,3);
t307 = -Icges(7,6) + Icges(6,6);
t306 = -Icges(6,3) - Icges(7,2);
t305 = rSges(7,1) + pkin(5);
t304 = rSges(7,3) + qJ(6);
t227 = pkin(10) + qJ(4);
t221 = cos(t227);
t228 = qJ(1) + pkin(9);
t222 = cos(t228);
t234 = cos(qJ(5));
t272 = t222 * t234;
t220 = sin(t228);
t232 = sin(qJ(5));
t275 = t220 * t232;
t171 = t221 * t275 + t272;
t273 = t222 * t232;
t274 = t220 * t234;
t172 = t221 * t274 - t273;
t219 = sin(t227);
t277 = t219 * t220;
t303 = t308 * t171 + t310 * t172 - t307 * t277;
t173 = t221 * t273 - t274;
t174 = t221 * t272 + t275;
t276 = t219 * t222;
t302 = t308 * t173 + t310 * t174 - t307 * t276;
t301 = -t307 * t171 + t309 * t172 - t306 * t277;
t300 = -t307 * t173 + t309 * t174 - t306 * t276;
t299 = t310 * t171 + t311 * t172 + t309 * t277;
t298 = t310 * t173 + t311 * t174 + t309 * t276;
t297 = t307 * t221 + (t308 * t232 + t310 * t234) * t219;
t296 = t306 * t221 + (-t307 * t232 + t309 * t234) * t219;
t295 = -t309 * t221 + (t310 * t232 + t311 * t234) * t219;
t233 = sin(qJ(1));
t235 = cos(qJ(1));
t294 = Icges(2,5) * t233 + Icges(3,5) * t220 + Icges(2,6) * t235 + Icges(3,6) * t222;
t293 = Icges(2,5) * t235 + Icges(3,5) * t222 - Icges(2,6) * t233 - Icges(3,6) * t220;
t288 = pkin(1) * t233;
t287 = pkin(1) * t235;
t229 = sin(pkin(10));
t286 = pkin(3) * t229;
t230 = cos(pkin(10));
t285 = pkin(3) * t230;
t284 = -pkin(6) - qJ(2);
t283 = Icges(2,4) * t233;
t282 = Icges(3,4) * t220;
t281 = Icges(4,4) * t229;
t280 = Icges(4,4) * t230;
t279 = Icges(5,4) * t219;
t278 = Icges(5,4) * t221;
t270 = rSges(7,2) * t277 + t304 * t171 + t305 * t172;
t269 = rSges(7,2) * t276 + t304 * t173 + t305 * t174;
t268 = -rSges(7,2) * t221 + (t304 * t232 + t305 * t234) * t219;
t267 = qJD(5) * t219;
t223 = V_base(6) + qJD(1);
t266 = t223 * t287 + V_base(2);
t265 = V_base(5) * pkin(6) + V_base(1);
t201 = qJD(4) * t220 + V_base(4);
t188 = pkin(2) * t220 - qJ(3) * t222;
t262 = -t188 - t288;
t190 = pkin(2) * t222 + qJ(3) * t220;
t261 = -t190 - t287;
t260 = V_base(5) * qJ(2) + t265;
t259 = V_base(4) * t288 + qJD(2) + V_base(3);
t136 = -pkin(7) * t222 + t220 * t285;
t258 = -t136 + t262;
t257 = qJD(3) * t220 + t260;
t256 = pkin(4) * t221 + pkin(8) * t219;
t255 = V_base(4) * t188 + t259;
t200 = -qJD(4) * t222 + V_base(5);
t254 = rSges(4,1) * t230 - rSges(4,2) * t229;
t253 = rSges(5,1) * t221 - rSges(5,2) * t219;
t252 = Icges(4,1) * t230 - t281;
t251 = Icges(5,1) * t221 - t279;
t250 = -Icges(4,2) * t229 + t280;
t249 = -Icges(5,2) * t219 + t278;
t248 = Icges(4,5) * t230 - Icges(4,6) * t229;
t247 = Icges(5,5) * t221 - Icges(5,6) * t219;
t246 = V_base(5) * t286 + t257;
t245 = -qJD(3) * t222 + t223 * t190 + t266;
t244 = (-Icges(5,3) * t222 + t220 * t247) * t200 + (Icges(5,3) * t220 + t222 * t247) * t201 + (Icges(5,5) * t219 + Icges(5,6) * t221) * t223;
t243 = (-Icges(4,3) * t222 + t220 * t248) * V_base(5) + (Icges(4,3) * t220 + t222 * t248) * V_base(4) + (Icges(4,5) * t229 + Icges(4,6) * t230) * t223;
t137 = pkin(7) * t220 + t222 * t285;
t242 = V_base(4) * t136 + (-t137 + t261) * V_base(5) + t255;
t167 = t256 * t220;
t192 = pkin(4) * t219 - pkin(8) * t221;
t241 = t200 * t192 + (-t167 + t258) * t223 + t246;
t240 = t223 * t137 + (t284 - t286) * V_base(4) + t245;
t168 = t256 * t222;
t239 = t201 * t167 - t200 * t168 + t242;
t238 = t223 * t168 - t192 * t201 + t240;
t140 = -Icges(5,6) * t222 + t220 * t249;
t141 = Icges(5,6) * t220 + t222 * t249;
t142 = -Icges(5,5) * t222 + t220 * t251;
t143 = Icges(5,5) * t220 + t222 * t251;
t181 = Icges(5,2) * t221 + t279;
t184 = Icges(5,1) * t219 + t278;
t237 = (-t141 * t219 + t143 * t221) * t201 + (-t140 * t219 + t142 * t221) * t200 + (-t181 * t219 + t184 * t221) * t223;
t150 = -Icges(4,6) * t222 + t220 * t250;
t151 = Icges(4,6) * t220 + t222 * t250;
t152 = -Icges(4,5) * t222 + t220 * t252;
t153 = Icges(4,5) * t220 + t222 * t252;
t198 = Icges(4,2) * t230 + t281;
t199 = Icges(4,1) * t229 + t280;
t236 = (-t151 * t229 + t153 * t230) * V_base(4) + (-t150 * t229 + t152 * t230) * V_base(5) + (-t198 * t229 + t199 * t230) * t223;
t225 = Icges(2,4) * t235;
t217 = Icges(3,4) * t222;
t210 = rSges(2,1) * t235 - t233 * rSges(2,2);
t209 = t233 * rSges(2,1) + rSges(2,2) * t235;
t208 = Icges(2,1) * t235 - t283;
t207 = Icges(2,1) * t233 + t225;
t206 = -Icges(2,2) * t233 + t225;
t205 = Icges(2,2) * t235 + t283;
t202 = rSges(4,1) * t229 + rSges(4,2) * t230;
t196 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t195 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t194 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t193 = -qJD(5) * t221 + t223;
t191 = rSges(3,1) * t222 - rSges(3,2) * t220;
t189 = rSges(3,1) * t220 + rSges(3,2) * t222;
t187 = rSges(5,1) * t219 + rSges(5,2) * t221;
t186 = Icges(3,1) * t222 - t282;
t185 = Icges(3,1) * t220 + t217;
t183 = -Icges(3,2) * t220 + t217;
t182 = Icges(3,2) * t222 + t282;
t170 = t222 * t267 + t201;
t169 = t220 * t267 + t200;
t165 = -rSges(6,3) * t221 + (rSges(6,1) * t234 - rSges(6,2) * t232) * t219;
t163 = V_base(5) * rSges(2,3) - t209 * t223 + t265;
t162 = t210 * t223 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = rSges(4,3) * t220 + t222 * t254;
t154 = -rSges(4,3) * t222 + t220 * t254;
t147 = t209 * V_base(4) - t210 * V_base(5) + V_base(3);
t145 = rSges(5,3) * t220 + t222 * t253;
t144 = -rSges(5,3) * t222 + t220 * t253;
t132 = V_base(5) * rSges(3,3) + (-t189 - t288) * t223 + t260;
t131 = t191 * t223 + (-rSges(3,3) + t284) * V_base(4) + t266;
t128 = V_base(4) * t189 + (-t191 - t287) * V_base(5) + t259;
t127 = rSges(6,1) * t174 - rSges(6,2) * t173 + rSges(6,3) * t276;
t125 = rSges(6,1) * t172 - rSges(6,2) * t171 + rSges(6,3) * t277;
t111 = t202 * V_base(5) + (-t154 + t262) * t223 + t257;
t110 = t155 * t223 + (-t202 + t284) * V_base(4) + t245;
t109 = V_base(4) * t154 + (-t155 + t261) * V_base(5) + t255;
t108 = t187 * t200 + (-t144 + t258) * t223 + t246;
t107 = t145 * t223 - t187 * t201 + t240;
t106 = t201 * t144 - t200 * t145 + t242;
t105 = -t125 * t193 + t165 * t169 + t241;
t104 = t127 * t193 - t165 * t170 + t238;
t103 = t170 * t125 - t169 * t127 + t239;
t102 = qJD(6) * t173 + t169 * t268 - t193 * t270 + t241;
t101 = qJD(6) * t171 - t170 * t268 + t193 * t269 + t238;
t100 = qJD(6) * t219 * t232 - t169 * t269 + t170 * t270 + t239;
t1 = t201 * (t244 * t220 + t237 * t222) / 0.2e1 + t200 * (t237 * t220 - t244 * t222) / 0.2e1 + m(1) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t128 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + ((t171 * t297 + t172 * t295 + t277 * t296) * t193 + (t171 * t302 + t172 * t298 + t277 * t300) * t170 + (t303 * t171 + t299 * t172 + t301 * t277) * t169) * t169 / 0.2e1 + ((t173 * t297 + t174 * t295 + t276 * t296) * t193 + (t302 * t173 + t298 * t174 + t300 * t276) * t170 + (t173 * t303 + t299 * t174 + t301 * t276) * t169) * t170 / 0.2e1 + ((-t169 * t301 - t170 * t300 - t193 * t296) * t221 + ((t232 * t297 + t234 * t295) * t193 + (t232 * t302 + t234 * t298) * t170 + (t232 * t303 + t299 * t234) * t169) * t219) * t193 / 0.2e1 + ((t141 * t221 + t143 * t219) * t201 + (t140 * t221 + t142 * t219) * t200 + (t150 * t230 + t152 * t229 + t294) * V_base(5) + (t151 * t230 + t153 * t229 + t293) * V_base(4) + (t181 * t221 + t184 * t219 + t198 * t230 + t199 * t229 + Icges(2,3) + Icges(3,3)) * t223) * t223 / 0.2e1 + (t243 * t220 + t236 * t222 + t293 * t223 + (-t182 * t220 + t185 * t222 - t233 * t205 + t207 * t235 + Icges(1,4)) * V_base(5) + (-t183 * t220 + t186 * t222 - t233 * t206 + t208 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t236 * t220 - t243 * t222 + t294 * t223 + (t182 * t222 + t185 * t220 + t205 * t235 + t233 * t207 + Icges(1,2)) * V_base(5) + (t183 * t222 + t186 * t220 + t206 * t235 + t233 * t208 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
