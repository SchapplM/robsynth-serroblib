% Calculate kinetic energy for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:42
% EndTime: 2019-03-09 05:58:45
% DurationCPUTime: 3.23s
% Computational Cost: add. (2270->322), mult. (2127->462), div. (0->0), fcn. (2019->10), ass. (0->158)
t321 = Icges(6,1) + Icges(7,1);
t320 = Icges(6,4) + Icges(7,4);
t319 = -Icges(7,5) - Icges(6,5);
t318 = Icges(6,2) + Icges(7,2);
t317 = -Icges(7,6) - Icges(6,6);
t316 = -Icges(7,3) - Icges(6,3);
t241 = qJ(1) + pkin(10);
t231 = sin(t241);
t232 = cos(t241);
t242 = qJ(4) + qJ(5);
t235 = cos(t242);
t234 = sin(t242);
t247 = cos(qJ(3));
t286 = t234 * t247;
t165 = -t231 * t286 - t232 * t235;
t285 = t235 * t247;
t166 = t231 * t285 - t232 * t234;
t244 = sin(qJ(3));
t289 = t231 * t244;
t315 = -t317 * t165 - t319 * t166 - t316 * t289;
t167 = t231 * t235 - t232 * t286;
t168 = t231 * t234 + t232 * t285;
t287 = t232 * t244;
t314 = -t317 * t167 - t319 * t168 - t316 * t287;
t313 = t318 * t165 + t320 * t166 - t317 * t289;
t312 = t318 * t167 + t320 * t168 - t317 * t287;
t311 = t320 * t165 + t321 * t166 - t319 * t289;
t310 = t320 * t167 + t321 * t168 - t319 * t287;
t309 = t316 * t247 + (t317 * t234 - t319 * t235) * t244;
t308 = t317 * t247 + (-t318 * t234 + t320 * t235) * t244;
t307 = t319 * t247 + (-t320 * t234 + t321 * t235) * t244;
t245 = sin(qJ(1));
t300 = pkin(1) * t245;
t248 = cos(qJ(1));
t299 = pkin(1) * t248;
t246 = cos(qJ(4));
t297 = t246 * pkin(4);
t296 = -pkin(6) - qJ(2);
t294 = Icges(2,4) * t245;
t293 = Icges(3,4) * t231;
t292 = Icges(4,4) * t244;
t291 = Icges(4,4) * t247;
t243 = sin(qJ(4));
t290 = t231 * t243;
t288 = t232 * t243;
t284 = t243 * t247;
t283 = t246 * t247;
t279 = pkin(5) * t235;
t259 = qJ(6) * t244 + t247 * t279;
t269 = pkin(5) * t234;
t282 = rSges(7,1) * t166 + rSges(7,2) * t165 + rSges(7,3) * t289 + t231 * t259 - t232 * t269;
t281 = rSges(7,1) * t168 + rSges(7,2) * t167 + rSges(7,3) * t287 + t231 * t269 + t232 * t259;
t280 = (-qJ(6) - rSges(7,3)) * t247 + (rSges(7,1) * t235 - rSges(7,2) * t234 + t279) * t244;
t277 = qJD(4) * t244;
t276 = qJD(5) * t244;
t275 = qJD(6) * t244;
t233 = V_base(6) + qJD(1);
t274 = t233 * t299 + V_base(2);
t273 = V_base(5) * pkin(6) + V_base(1);
t209 = qJD(3) * t231 + V_base(4);
t201 = pkin(2) * t231 - pkin(7) * t232;
t270 = -t201 - t300;
t268 = V_base(5) * qJ(2) + t273;
t267 = V_base(4) * t300 + qJD(2) + V_base(3);
t179 = t232 * t277 + t209;
t266 = pkin(3) * t247 + pkin(8) * t244;
t208 = -qJD(3) * t232 + V_base(5);
t265 = rSges(4,1) * t247 - rSges(4,2) * t244;
t264 = Icges(4,1) * t247 - t292;
t263 = -Icges(4,2) * t244 + t291;
t262 = Icges(4,5) * t247 - Icges(4,6) * t244;
t178 = t231 * t277 + t208;
t261 = pkin(9) * t244 + t247 * t297;
t260 = (-Icges(4,3) * t232 + t231 * t262) * t208 + (Icges(4,3) * t231 + t232 * t262) * t209 + (Icges(4,5) * t244 + Icges(4,6) * t247) * t233;
t202 = pkin(2) * t232 + pkin(7) * t231;
t258 = t233 * t202 + t296 * V_base(4) + t274;
t189 = t266 * t231;
t225 = t244 * pkin(3) - t247 * pkin(8);
t257 = t208 * t225 + (-t189 + t270) * t233 + t268;
t256 = V_base(4) * t201 + (-t202 - t299) * V_base(5) + t267;
t190 = t266 * t232;
t255 = t233 * t190 - t209 * t225 + t258;
t142 = -pkin(4) * t288 + t231 * t261;
t164 = -pkin(9) * t247 + t244 * t297;
t221 = -qJD(4) * t247 + t233;
t254 = -t142 * t221 + t178 * t164 + t257;
t253 = t209 * t189 - t190 * t208 + t256;
t143 = pkin(4) * t290 + t232 * t261;
t252 = t221 * t143 - t164 * t179 + t255;
t251 = t179 * t142 - t143 * t178 + t253;
t155 = -Icges(4,6) * t232 + t231 * t263;
t156 = Icges(4,6) * t231 + t232 * t263;
t157 = -Icges(4,5) * t232 + t231 * t264;
t158 = Icges(4,5) * t231 + t232 * t264;
t213 = Icges(4,2) * t247 + t292;
t216 = Icges(4,1) * t244 + t291;
t250 = (-t156 * t244 + t158 * t247) * t209 + (-t155 * t244 + t157 * t247) * t208 + (-t213 * t244 + t216 * t247) * t233;
t237 = Icges(2,4) * t248;
t229 = Icges(3,4) * t232;
t224 = rSges(2,1) * t248 - rSges(2,2) * t245;
t223 = rSges(2,1) * t245 + rSges(2,2) * t248;
t222 = rSges(4,1) * t244 + rSges(4,2) * t247;
t218 = Icges(2,1) * t248 - t294;
t217 = Icges(2,1) * t245 + t237;
t215 = -Icges(2,2) * t245 + t237;
t214 = Icges(2,2) * t248 + t294;
t206 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t205 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t204 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t200 = rSges(3,1) * t232 - rSges(3,2) * t231;
t199 = rSges(3,1) * t231 + rSges(3,2) * t232;
t198 = (-qJD(4) - qJD(5)) * t247 + t233;
t197 = Icges(3,1) * t232 - t293;
t196 = Icges(3,1) * t231 + t229;
t195 = -Icges(3,2) * t231 + t229;
t194 = Icges(3,2) * t232 + t293;
t187 = -rSges(5,3) * t247 + (rSges(5,1) * t246 - rSges(5,2) * t243) * t244;
t186 = -Icges(5,5) * t247 + (Icges(5,1) * t246 - Icges(5,4) * t243) * t244;
t185 = -Icges(5,6) * t247 + (Icges(5,4) * t246 - Icges(5,2) * t243) * t244;
t184 = -Icges(5,3) * t247 + (Icges(5,5) * t246 - Icges(5,6) * t243) * t244;
t183 = t232 * t283 + t290;
t182 = t231 * t246 - t232 * t284;
t181 = t231 * t283 - t288;
t180 = -t231 * t284 - t232 * t246;
t176 = -rSges(6,3) * t247 + (rSges(6,1) * t235 - rSges(6,2) * t234) * t244;
t162 = rSges(4,3) * t231 + t232 * t265;
t161 = -rSges(4,3) * t232 + t231 * t265;
t160 = V_base(5) * rSges(2,3) - t223 * t233 + t273;
t159 = t224 * t233 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t152 = t223 * V_base(4) - t224 * V_base(5) + V_base(3);
t151 = t232 * t276 + t179;
t150 = t231 * t276 + t178;
t147 = V_base(5) * rSges(3,3) + (-t199 - t300) * t233 + t268;
t146 = t200 * t233 + (-rSges(3,3) + t296) * V_base(4) + t274;
t144 = t199 * V_base(4) + (-t200 - t299) * V_base(5) + t267;
t141 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t287;
t140 = rSges(5,1) * t181 + rSges(5,2) * t180 + rSges(5,3) * t289;
t139 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t287;
t138 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t289;
t137 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t287;
t136 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t289;
t135 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t287;
t134 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t289;
t133 = rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t287;
t131 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t289;
t113 = t208 * t222 + (-t161 + t270) * t233 + t268;
t112 = t162 * t233 - t209 * t222 + t258;
t111 = t161 * t209 - t162 * t208 + t256;
t110 = -t140 * t221 + t178 * t187 + t257;
t109 = t141 * t221 - t179 * t187 + t255;
t108 = t140 * t179 - t141 * t178 + t253;
t107 = -t131 * t198 + t150 * t176 + t254;
t106 = t133 * t198 - t151 * t176 + t252;
t105 = t131 * t151 - t133 * t150 + t251;
t104 = t150 * t280 - t198 * t282 + t232 * t275 + t254;
t103 = -t151 * t280 + t198 * t281 + t231 * t275 + t252;
t102 = -qJD(6) * t247 - t150 * t281 + t151 * t282 + t251;
t1 = t179 * ((t135 * t287 + t182 * t137 + t183 * t139) * t179 + (t134 * t287 + t136 * t182 + t138 * t183) * t178 + (t182 * t185 + t183 * t186 + t184 * t287) * t221) / 0.2e1 + t178 * ((t135 * t289 + t137 * t180 + t139 * t181) * t179 + (t134 * t289 + t180 * t136 + t181 * t138) * t178 + (t180 * t185 + t181 * t186 + t184 * t289) * t221) / 0.2e1 + t209 * (t260 * t231 + t250 * t232) / 0.2e1 + t208 * (t250 * t231 - t260 * t232) / 0.2e1 + m(7) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(6) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + t221 * ((-t134 * t178 - t135 * t179 - t184 * t221) * t247 + ((-t137 * t243 + t139 * t246) * t179 + (-t136 * t243 + t138 * t246) * t178 + (-t185 * t243 + t186 * t246) * t221) * t244) / 0.2e1 + m(3) * (t144 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(2) * (t152 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(1) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + ((t165 * t308 + t166 * t307 + t289 * t309) * t198 + (t165 * t312 + t166 * t310 + t289 * t314) * t151 + (t313 * t165 + t311 * t166 + t315 * t289) * t150) * t150 / 0.2e1 + ((t167 * t308 + t168 * t307 + t287 * t309) * t198 + (t312 * t167 + t310 * t168 + t314 * t287) * t151 + (t313 * t167 + t311 * t168 + t287 * t315) * t150) * t151 / 0.2e1 + ((-t150 * t315 - t314 * t151 - t309 * t198) * t247 + ((-t234 * t308 + t235 * t307) * t198 + (-t234 * t312 + t235 * t310) * t151 + (-t234 * t313 + t235 * t311) * t150) * t244) * t198 / 0.2e1 + ((t156 * t247 + t158 * t244) * t209 + (t155 * t247 + t157 * t244) * t208 + (t213 * t247 + t216 * t244 + Icges(2,3) + Icges(3,3)) * t233) * t233 / 0.2e1 + ((-t194 * t231 + t196 * t232 - t214 * t245 + t217 * t248 + Icges(1,4)) * V_base(5) + (-t195 * t231 + t197 * t232 - t215 * t245 + t218 * t248 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t194 * t232 + t196 * t231 + t214 * t248 + t217 * t245 + Icges(1,2)) * V_base(5) + (t195 * t232 + t197 * t231 + t215 * t248 + t218 * t245 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t233 * (Icges(2,5) * t245 + Icges(3,5) * t231 + Icges(2,6) * t248 + Icges(3,6) * t232) + V_base(4) * t233 * (Icges(2,5) * t248 + Icges(3,5) * t232 - Icges(2,6) * t245 - Icges(3,6) * t231) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
