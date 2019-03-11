% Calculate kinetic energy for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:52
% EndTime: 2019-03-09 03:47:56
% DurationCPUTime: 4.30s
% Computational Cost: add. (1887->333), mult. (2303->471), div. (0->0), fcn. (2319->10), ass. (0->165)
t332 = Icges(4,4) - Icges(5,5);
t331 = Icges(4,1) + Icges(5,1);
t330 = Icges(4,2) + Icges(5,3);
t249 = pkin(10) + qJ(3);
t241 = sin(t249);
t329 = t332 * t241;
t242 = cos(t249);
t328 = t332 * t242;
t327 = Icges(5,4) + Icges(4,5);
t326 = Icges(4,6) - Icges(5,6);
t325 = t330 * t241 - t328;
t324 = t331 * t242 - t329;
t323 = Icges(5,2) + Icges(4,3);
t255 = sin(qJ(1));
t257 = cos(qJ(1));
t322 = t325 * t255 + t326 * t257;
t321 = -t326 * t255 + t325 * t257;
t320 = t324 * t255 - t327 * t257;
t319 = t327 * t255 + t324 * t257;
t318 = -t330 * t242 - t329;
t317 = t331 * t241 + t328;
t316 = -t326 * t241 + t327 * t242;
t237 = -qJD(3) * t257 + V_base(5);
t238 = qJD(3) * t255 + V_base(4);
t243 = V_base(6) + qJD(1);
t315 = (t318 * t241 + t317 * t242) * t243 + (t321 * t241 + t319 * t242) * t238 + (t322 * t241 + t320 * t242) * t237;
t314 = (t327 * t241 + t326 * t242) * t243 + (t323 * t255 + t316 * t257) * t238 + (t316 * t255 - t323 * t257) * t237;
t310 = cos(qJ(5));
t250 = sin(pkin(10));
t309 = pkin(2) * t250;
t308 = pkin(4) * t241;
t251 = cos(pkin(10));
t307 = pkin(2) * t251;
t306 = Icges(2,4) * t255;
t305 = Icges(3,4) * t250;
t304 = Icges(3,4) * t251;
t299 = t242 * t255;
t298 = t242 * t257;
t161 = -pkin(7) * t257 + t255 * t307;
t233 = t255 * pkin(1) - qJ(2) * t257;
t296 = -t161 - t233;
t295 = qJD(4) * t241;
t294 = V_base(4) * t233 + V_base(3);
t293 = V_base(5) * pkin(6) + V_base(1);
t282 = pkin(3) * t242 + qJ(4) * t241;
t196 = t282 * t255;
t290 = -t196 + t296;
t289 = t241 * t310;
t288 = qJD(2) * t255 + t293;
t208 = pkin(4) * t299 + pkin(8) * t257;
t287 = -t208 + t290;
t286 = V_base(5) * t309 + t288;
t285 = rSges(3,1) * t251 - rSges(3,2) * t250;
t284 = rSges(4,1) * t242 - rSges(4,2) * t241;
t283 = rSges(5,1) * t242 + rSges(5,3) * t241;
t281 = Icges(3,1) * t251 - t305;
t278 = -Icges(3,2) * t250 + t304;
t275 = Icges(3,5) * t251 - Icges(3,6) * t250;
t235 = pkin(1) * t257 + t255 * qJ(2);
t272 = -qJD(2) * t257 + t243 * t235 + V_base(2);
t213 = qJD(5) * t257 + t237;
t214 = -qJD(5) * t255 + t238;
t254 = sin(qJ(5));
t198 = t241 * t254 + t242 * t310;
t210 = pkin(3) * t241 - qJ(4) * t242;
t271 = t237 * t210 + t257 * t295 + t286;
t270 = t237 * t308 + t271;
t162 = pkin(7) * t255 + t257 * t307;
t267 = V_base(4) * t161 + (-t162 - t235) * V_base(5) + t294;
t266 = (-Icges(3,3) * t257 + t255 * t275) * V_base(5) + (Icges(3,3) * t255 + t257 * t275) * V_base(4) + (Icges(3,5) * t250 + Icges(3,6) * t251) * t243;
t265 = t243 * t162 + (-pkin(6) - t309) * V_base(4) + t272;
t264 = -qJD(4) * t242 + t238 * t196 + t267;
t197 = t282 * t257;
t263 = t243 * t197 + t255 * t295 + t265;
t209 = pkin(4) * t298 - t255 * pkin(8);
t262 = t238 * t208 + (-t197 - t209) * t237 + t264;
t261 = t243 * t209 + (-t210 - t308) * t238 + t263;
t189 = -Icges(3,6) * t257 + t255 * t278;
t190 = Icges(3,6) * t255 + t257 * t278;
t191 = -Icges(3,5) * t257 + t255 * t281;
t192 = Icges(3,5) * t255 + t257 * t281;
t222 = Icges(3,2) * t251 + t305;
t223 = Icges(3,1) * t250 + t304;
t258 = (-t190 * t250 + t192 * t251) * V_base(4) + (-t189 * t250 + t191 * t251) * V_base(5) + (-t222 * t250 + t223 * t251) * t243;
t256 = cos(qJ(6));
t253 = sin(qJ(6));
t247 = Icges(2,4) * t257;
t236 = rSges(2,1) * t257 - t255 * rSges(2,2);
t234 = t255 * rSges(2,1) + rSges(2,2) * t257;
t230 = Icges(2,1) * t257 - t306;
t229 = Icges(2,1) * t255 + t247;
t228 = -Icges(2,2) * t255 + t247;
t227 = Icges(2,2) * t257 + t306;
t226 = Icges(2,5) * t257 - Icges(2,6) * t255;
t225 = Icges(2,5) * t255 + Icges(2,6) * t257;
t224 = rSges(3,1) * t250 + rSges(3,2) * t251;
t218 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t217 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t216 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t212 = rSges(4,1) * t241 + rSges(4,2) * t242;
t211 = rSges(5,1) * t241 - rSges(5,3) * t242;
t199 = -t242 * t254 + t289;
t194 = t255 * rSges(3,3) + t257 * t285;
t193 = -rSges(3,3) * t257 + t255 * t285;
t186 = t198 * t257;
t185 = t254 * t298 - t257 * t289;
t184 = t198 * t255;
t183 = t254 * t299 - t255 * t289;
t180 = t255 * rSges(4,3) + t257 * t284;
t179 = t255 * rSges(5,2) + t257 * t283;
t178 = -rSges(4,3) * t257 + t255 * t284;
t177 = -rSges(5,2) * t257 + t255 * t283;
t176 = qJD(6) * t198 + t243;
t160 = V_base(5) * rSges(2,3) - t234 * t243 + t293;
t159 = t236 * t243 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t157 = t234 * V_base(4) - t236 * V_base(5) + V_base(3);
t154 = t186 * t256 - t253 * t255;
t153 = -t186 * t253 - t255 * t256;
t152 = t184 * t256 + t253 * t257;
t151 = -t184 * t253 + t256 * t257;
t150 = qJD(6) * t185 + t214;
t149 = qJD(6) * t183 + t213;
t148 = pkin(5) * t199 + pkin(9) * t198;
t147 = rSges(6,1) * t199 - rSges(6,2) * t198;
t146 = Icges(6,1) * t199 - Icges(6,4) * t198;
t145 = Icges(6,4) * t199 - Icges(6,2) * t198;
t144 = Icges(6,5) * t199 - Icges(6,6) * t198;
t143 = pkin(5) * t186 + pkin(9) * t185;
t142 = pkin(5) * t184 + pkin(9) * t183;
t141 = rSges(6,1) * t186 - rSges(6,2) * t185 - rSges(6,3) * t255;
t140 = t184 * rSges(6,1) - t183 * rSges(6,2) + rSges(6,3) * t257;
t139 = Icges(6,1) * t186 - Icges(6,4) * t185 - Icges(6,5) * t255;
t138 = Icges(6,1) * t184 - Icges(6,4) * t183 + Icges(6,5) * t257;
t137 = Icges(6,4) * t186 - Icges(6,2) * t185 - Icges(6,6) * t255;
t136 = Icges(6,4) * t184 - Icges(6,2) * t183 + Icges(6,6) * t257;
t135 = Icges(6,5) * t186 - Icges(6,6) * t185 - Icges(6,3) * t255;
t134 = Icges(6,5) * t184 - Icges(6,6) * t183 + Icges(6,3) * t257;
t133 = rSges(7,3) * t198 + (rSges(7,1) * t256 - rSges(7,2) * t253) * t199;
t132 = Icges(7,5) * t198 + (Icges(7,1) * t256 - Icges(7,4) * t253) * t199;
t131 = Icges(7,6) * t198 + (Icges(7,4) * t256 - Icges(7,2) * t253) * t199;
t130 = Icges(7,3) * t198 + (Icges(7,5) * t256 - Icges(7,6) * t253) * t199;
t129 = t224 * V_base(5) + (-t193 - t233) * t243 + t288;
t128 = t243 * t194 + (-pkin(6) - t224) * V_base(4) + t272;
t127 = t193 * V_base(4) + (-t194 - t235) * V_base(5) + t294;
t126 = rSges(7,1) * t154 + rSges(7,2) * t153 + rSges(7,3) * t185;
t125 = rSges(7,1) * t152 + rSges(7,2) * t151 + rSges(7,3) * t183;
t124 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t185;
t123 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t183;
t122 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t185;
t121 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t183;
t120 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t185;
t119 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t183;
t118 = t212 * t237 + (-t178 + t296) * t243 + t286;
t117 = t243 * t180 - t238 * t212 + t265;
t116 = t178 * t238 - t180 * t237 + t267;
t115 = t211 * t237 + (-t177 + t290) * t243 + t271;
t114 = t243 * t179 + (-t210 - t211) * t238 + t263;
t113 = t177 * t238 + (-t179 - t197) * t237 + t264;
t112 = t147 * t213 + (-t140 + t287) * t243 + t270;
t111 = t243 * t141 - t214 * t147 + t261;
t110 = t140 * t214 - t141 * t213 + t262;
t109 = -t125 * t176 + t133 * t149 + t148 * t213 + (-t142 + t287) * t243 + t270;
t108 = t176 * t126 - t150 * t133 + t243 * t143 - t214 * t148 + t261;
t107 = t125 * t150 - t126 * t149 + t142 * t214 - t143 * t213 + t262;
t1 = t213 * ((t135 * t257 - t183 * t137 + t184 * t139) * t214 + (t257 * t134 - t183 * t136 + t184 * t138) * t213 + (t144 * t257 - t183 * t145 + t184 * t146) * t243) / 0.2e1 + t214 * ((-t255 * t135 - t185 * t137 + t186 * t139) * t214 + (-t134 * t255 - t136 * t185 + t138 * t186) * t213 + (-t144 * t255 - t145 * t185 + t146 * t186) * t243) / 0.2e1 + m(1) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + t149 * ((t120 * t183 + t122 * t151 + t124 * t152) * t150 + (t183 * t119 + t151 * t121 + t152 * t123) * t149 + (t130 * t183 + t131 * t151 + t132 * t152) * t176) / 0.2e1 + t150 * ((t185 * t120 + t153 * t122 + t154 * t124) * t150 + (t119 * t185 + t121 * t153 + t123 * t154) * t149 + (t130 * t185 + t131 * t153 + t132 * t154) * t176) / 0.2e1 + m(2) * (t157 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t176 * ((t119 * t149 + t120 * t150 + t130 * t176) * t198 + ((-t122 * t253 + t124 * t256) * t150 + (-t121 * t253 + t123 * t256) * t149 + (-t131 * t253 + t132 * t256) * t176) * t199) / 0.2e1 + (t315 * t255 - t314 * t257) * t237 / 0.2e1 + (t314 * t255 + t315 * t257) * t238 / 0.2e1 + (t226 * t243 + t266 * t255 + t258 * t257 + (-t255 * t227 + t229 * t257 + Icges(1,4)) * V_base(5) + (-t255 * t228 + t230 * t257 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t225 * t243 + t258 * t255 - t266 * t257 + (t227 * t257 + t255 * t229 + Icges(1,2)) * V_base(5) + (t228 * t257 + t255 * t230 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t137 * t198 + t139 * t199) * t214 + (-t136 * t198 + t138 * t199) * t213 + (t189 * t251 + t191 * t250 + t225) * V_base(5) + (t190 * t251 + t192 * t250 + t226) * V_base(4) + (t319 * t241 - t321 * t242) * t238 + (t320 * t241 - t322 * t242) * t237 + (-t145 * t198 + t146 * t199 + t222 * t251 + t223 * t250 + t317 * t241 - t318 * t242 + Icges(2,3)) * t243) * t243 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
