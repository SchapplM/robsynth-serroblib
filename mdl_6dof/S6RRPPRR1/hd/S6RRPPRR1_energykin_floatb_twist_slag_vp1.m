% Calculate kinetic energy for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:43
% EndTime: 2019-03-09 08:45:46
% DurationCPUTime: 3.99s
% Computational Cost: add. (1913->326), mult. (2329->457), div. (0->0), fcn. (2345->10), ass. (0->161)
t333 = Icges(4,4) - Icges(5,5);
t332 = Icges(4,1) + Icges(5,1);
t331 = Icges(4,2) + Icges(5,3);
t249 = qJ(2) + pkin(10);
t241 = sin(t249);
t330 = t333 * t241;
t242 = cos(t249);
t329 = t333 * t242;
t328 = Icges(5,4) + Icges(4,5);
t327 = Icges(4,6) - Icges(5,6);
t326 = t331 * t241 - t329;
t325 = t332 * t242 - t330;
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t324 = t254 * t326 + t257 * t327;
t323 = -t254 * t327 + t257 * t326;
t322 = t325 * t254 - t257 * t328;
t321 = t254 * t328 + t325 * t257;
t320 = -t331 * t242 - t330;
t319 = t332 * t241 + t329;
t318 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t253 = sin(qJ(2));
t256 = cos(qJ(2));
t317 = Icges(3,5) * t256 - Icges(3,6) * t253 - t327 * t241 + t242 * t328;
t304 = Icges(3,4) * t256;
t279 = -Icges(3,2) * t253 + t304;
t189 = -Icges(3,6) * t257 + t254 * t279;
t190 = Icges(3,6) * t254 + t257 * t279;
t305 = Icges(3,4) * t253;
t282 = Icges(3,1) * t256 - t305;
t191 = -Icges(3,5) * t257 + t254 * t282;
t192 = Icges(3,5) * t254 + t257 * t282;
t225 = Icges(3,2) * t256 + t305;
t228 = Icges(3,1) * t253 + t304;
t238 = -qJD(2) * t257 + V_base(5);
t239 = qJD(2) * t254 + V_base(4);
t243 = V_base(6) + qJD(1);
t316 = (-t225 * t253 + t228 * t256 + t241 * t320 + t242 * t319) * t243 + (-t190 * t253 + t192 * t256 + t241 * t323 + t242 * t321) * t239 + (-t189 * t253 + t191 * t256 + t241 * t324 + t322 * t242) * t238;
t315 = (Icges(3,5) * t253 + Icges(3,6) * t256 + t241 * t328 + t327 * t242) * t243 + (t254 * t318 + t257 * t317) * t239 + (t254 * t317 - t257 * t318) * t238;
t311 = cos(qJ(5));
t310 = pkin(2) * t253;
t309 = pkin(4) * t241;
t308 = pkin(2) * t256;
t306 = Icges(2,4) * t254;
t299 = t242 * t254;
t298 = t242 * t257;
t162 = -qJ(3) * t257 + t254 * t308;
t236 = t254 * pkin(1) - pkin(7) * t257;
t297 = -t162 - t236;
t163 = qJ(3) * t254 + t257 * t308;
t283 = pkin(3) * t242 + qJ(4) * t241;
t195 = t283 * t257;
t296 = -t163 - t195;
t295 = qJD(4) * t241;
t294 = V_base(5) * pkin(6) + V_base(1);
t194 = t283 * t254;
t291 = -t194 + t297;
t290 = t241 * t311;
t210 = pkin(3) * t241 - qJ(4) * t242;
t289 = -t210 - t310;
t208 = pkin(4) * t299 + pkin(8) * t257;
t288 = -t208 + t291;
t287 = qJD(3) * t254 + t238 * t310 + t294;
t286 = rSges(3,1) * t256 - rSges(3,2) * t253;
t285 = rSges(4,1) * t242 - rSges(4,2) * t241;
t284 = rSges(5,1) * t242 + rSges(5,3) * t241;
t213 = qJD(5) * t257 + t238;
t214 = -qJD(5) * t254 + t239;
t237 = pkin(1) * t257 + t254 * pkin(7);
t273 = -V_base(4) * pkin(6) + t243 * t237 + V_base(2);
t252 = sin(qJ(5));
t198 = t241 * t252 + t242 * t311;
t272 = V_base(4) * t236 - t237 * V_base(5) + V_base(3);
t271 = t238 * t210 + t257 * t295 + t287;
t270 = t239 * t162 + t272;
t269 = t238 * t309 + t271;
t265 = -qJD(3) * t257 + t243 * t163 + t273;
t264 = -qJD(4) * t242 + t239 * t194 + t270;
t263 = t243 * t195 + t254 * t295 + t265;
t209 = pkin(4) * t298 - t254 * pkin(8);
t262 = t239 * t208 + (-t209 + t296) * t238 + t264;
t261 = t243 * t209 + (t289 - t309) * t239 + t263;
t255 = cos(qJ(6));
t251 = sin(qJ(6));
t247 = Icges(2,4) * t257;
t235 = rSges(2,1) * t257 - t254 * rSges(2,2);
t234 = t254 * rSges(2,1) + rSges(2,2) * t257;
t233 = rSges(3,1) * t253 + rSges(3,2) * t256;
t230 = Icges(2,1) * t257 - t306;
t229 = Icges(2,1) * t254 + t247;
t227 = -Icges(2,2) * t254 + t247;
t226 = Icges(2,2) * t257 + t306;
t219 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t218 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t217 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t212 = rSges(4,1) * t241 + rSges(4,2) * t242;
t211 = rSges(5,1) * t241 - rSges(5,3) * t242;
t199 = -t242 * t252 + t290;
t197 = t254 * rSges(3,3) + t257 * t286;
t196 = -rSges(3,3) * t257 + t254 * t286;
t186 = t198 * t257;
t185 = t252 * t298 - t257 * t290;
t184 = t198 * t254;
t183 = t252 * t299 - t254 * t290;
t180 = t254 * rSges(4,3) + t257 * t285;
t179 = t254 * rSges(5,2) + t257 * t284;
t178 = -rSges(4,3) * t257 + t254 * t285;
t177 = -rSges(5,2) * t257 + t254 * t284;
t176 = qJD(6) * t198 + t243;
t160 = V_base(5) * rSges(2,3) - t234 * t243 + t294;
t159 = t235 * t243 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t157 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t155 = t186 * t255 - t251 * t254;
t154 = -t186 * t251 - t254 * t255;
t153 = t184 * t255 + t251 * t257;
t152 = -t184 * t251 + t255 * t257;
t150 = qJD(6) * t185 + t214;
t149 = qJD(6) * t183 + t213;
t148 = pkin(5) * t199 + pkin(9) * t198;
t147 = rSges(6,1) * t199 - rSges(6,2) * t198;
t146 = Icges(6,1) * t199 - Icges(6,4) * t198;
t145 = Icges(6,4) * t199 - Icges(6,2) * t198;
t144 = Icges(6,5) * t199 - Icges(6,6) * t198;
t143 = pkin(5) * t186 + pkin(9) * t185;
t142 = pkin(5) * t184 + pkin(9) * t183;
t141 = rSges(6,1) * t186 - rSges(6,2) * t185 - rSges(6,3) * t254;
t140 = t184 * rSges(6,1) - t183 * rSges(6,2) + rSges(6,3) * t257;
t139 = Icges(6,1) * t186 - Icges(6,4) * t185 - Icges(6,5) * t254;
t138 = Icges(6,1) * t184 - Icges(6,4) * t183 + Icges(6,5) * t257;
t137 = Icges(6,4) * t186 - Icges(6,2) * t185 - Icges(6,6) * t254;
t136 = Icges(6,4) * t184 - Icges(6,2) * t183 + Icges(6,6) * t257;
t135 = Icges(6,5) * t186 - Icges(6,6) * t185 - Icges(6,3) * t254;
t134 = Icges(6,5) * t184 - Icges(6,6) * t183 + Icges(6,3) * t257;
t133 = rSges(7,3) * t198 + (rSges(7,1) * t255 - rSges(7,2) * t251) * t199;
t132 = Icges(7,5) * t198 + (Icges(7,1) * t255 - Icges(7,4) * t251) * t199;
t131 = Icges(7,6) * t198 + (Icges(7,4) * t255 - Icges(7,2) * t251) * t199;
t130 = Icges(7,3) * t198 + (Icges(7,5) * t255 - Icges(7,6) * t251) * t199;
t129 = t233 * t238 + (-t196 - t236) * t243 + t294;
t128 = t197 * t243 - t233 * t239 + t273;
t127 = t196 * t239 - t197 * t238 + t272;
t126 = rSges(7,1) * t155 + rSges(7,2) * t154 + rSges(7,3) * t185;
t125 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t183;
t124 = Icges(7,1) * t155 + Icges(7,4) * t154 + Icges(7,5) * t185;
t123 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t183;
t122 = Icges(7,4) * t155 + Icges(7,2) * t154 + Icges(7,6) * t185;
t121 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t183;
t120 = Icges(7,5) * t155 + Icges(7,6) * t154 + Icges(7,3) * t185;
t119 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t183;
t118 = t212 * t238 + (-t178 + t297) * t243 + t287;
t117 = t243 * t180 + (-t212 - t310) * t239 + t265;
t116 = t178 * t239 + (-t163 - t180) * t238 + t270;
t115 = t211 * t238 + (-t177 + t291) * t243 + t271;
t114 = t243 * t179 + (-t211 + t289) * t239 + t263;
t113 = t177 * t239 + (-t179 + t296) * t238 + t264;
t112 = t147 * t213 + (-t140 + t288) * t243 + t269;
t111 = t243 * t141 - t214 * t147 + t261;
t110 = t140 * t214 - t141 * t213 + t262;
t109 = -t125 * t176 + t133 * t149 + t148 * t213 + (-t142 + t288) * t243 + t269;
t108 = t176 * t126 - t150 * t133 + t243 * t143 - t214 * t148 + t261;
t107 = t125 * t150 - t126 * t149 + t142 * t214 - t143 * t213 + t262;
t1 = m(1) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + t149 * ((t120 * t183 + t122 * t152 + t124 * t153) * t150 + (t183 * t119 + t152 * t121 + t153 * t123) * t149 + (t130 * t183 + t131 * t152 + t132 * t153) * t176) / 0.2e1 + t150 * ((t185 * t120 + t154 * t122 + t155 * t124) * t150 + (t119 * t185 + t121 * t154 + t123 * t155) * t149 + (t130 * t185 + t131 * t154 + t132 * t155) * t176) / 0.2e1 + m(2) * (t157 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t214 * ((-t254 * t135 - t185 * t137 + t186 * t139) * t214 + (-t134 * t254 - t136 * t185 + t138 * t186) * t213 + (-t144 * t254 - t145 * t185 + t146 * t186) * t243) / 0.2e1 + t213 * ((t135 * t257 - t183 * t137 + t184 * t139) * t214 + (t257 * t134 - t183 * t136 + t184 * t138) * t213 + (t144 * t257 - t183 * t145 + t184 * t146) * t243) / 0.2e1 + t176 * ((t119 * t149 + t120 * t150 + t130 * t176) * t198 + ((-t122 * t251 + t124 * t255) * t150 + (-t121 * t251 + t123 * t255) * t149 + (-t131 * t251 + t132 * t255) * t176) * t199) / 0.2e1 + ((-t254 * t226 + t229 * t257 + Icges(1,4)) * V_base(5) + (-t254 * t227 + t257 * t230 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t257 * t226 + t254 * t229 + Icges(1,2)) * V_base(5) + (t227 * t257 + t254 * t230 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t316 * t254 - t315 * t257) * t238 / 0.2e1 + (t315 * t254 + t316 * t257) * t239 / 0.2e1 + ((-t137 * t198 + t139 * t199) * t214 + (-t136 * t198 + t138 * t199) * t213 + (t190 * t256 + t192 * t253 + t241 * t321 - t242 * t323) * t239 + (t189 * t256 + t191 * t253 + t322 * t241 - t242 * t324) * t238 + (-t198 * t145 + t199 * t146 + t256 * t225 + t253 * t228 + t241 * t319 - t320 * t242 + Icges(2,3)) * t243) * t243 / 0.2e1 + t243 * V_base(4) * (Icges(2,5) * t257 - Icges(2,6) * t254) + t243 * V_base(5) * (Icges(2,5) * t254 + Icges(2,6) * t257) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
