% Calculate kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:00
% EndTime: 2022-01-23 08:59:02
% DurationCPUTime: 2.62s
% Computational Cost: add. (1469->362), mult. (2999->473), div. (0->0), fcn. (3389->10), ass. (0->174)
t237 = sin(pkin(8));
t240 = cos(pkin(7));
t244 = cos(qJ(1));
t239 = cos(pkin(8));
t242 = sin(qJ(1));
t279 = t242 * t239;
t195 = -t237 * t244 + t240 * t279;
t236 = sin(pkin(9));
t238 = sin(pkin(7));
t291 = cos(pkin(9));
t264 = t238 * t291;
t160 = t195 * t236 - t242 * t264;
t283 = t238 * t242;
t161 = t195 * t291 + t236 * t283;
t280 = t242 * t237;
t194 = t239 * t244 + t240 * t280;
t122 = Icges(5,5) * t161 - Icges(5,6) * t160 + Icges(5,3) * t194;
t148 = Icges(4,4) * t195 - Icges(4,2) * t194 + Icges(4,6) * t283;
t299 = t122 - t148;
t281 = t240 * t244;
t197 = t239 * t281 + t280;
t262 = t244 * t291;
t162 = t197 * t236 - t238 * t262;
t282 = t238 * t244;
t163 = t197 * t291 + t236 * t282;
t196 = t237 * t281 - t279;
t123 = Icges(5,5) * t163 - Icges(5,6) * t162 + Icges(5,3) * t196;
t149 = Icges(4,4) * t197 - Icges(4,2) * t196 + Icges(4,6) * t282;
t298 = t123 - t149;
t263 = t240 * t291;
t284 = t238 * t236;
t190 = t239 * t284 + t263;
t191 = -t240 * t236 + t239 * t264;
t287 = t237 * t238;
t142 = Icges(5,5) * t191 - Icges(5,6) * t190 + Icges(5,3) * t287;
t172 = -Icges(4,6) * t240 + (Icges(4,4) * t239 - Icges(4,2) * t237) * t238;
t297 = t142 - t172;
t294 = t238 * qJ(3) + pkin(1);
t205 = pkin(2) * t240 + t294;
t295 = -pkin(1) + t205;
t293 = t237 * qJ(4) + pkin(2);
t211 = pkin(2) * t238 - t240 * qJ(3);
t292 = -pkin(5) - t211;
t290 = Icges(2,4) * t242;
t289 = Icges(3,4) * t238;
t288 = Icges(3,4) * t240;
t241 = sin(qJ(5));
t286 = t237 * t241;
t243 = cos(qJ(5));
t285 = t237 * t243;
t260 = qJ(4) * t239 - qJ(2);
t204 = -t237 * pkin(3) + t260;
t278 = qJ(2) + t204;
t206 = pkin(3) * t239 + t293;
t174 = t206 * t240 + t294;
t207 = pkin(4) * t291 + t236 * pkin(6) + pkin(3);
t176 = t207 * t239 + t293;
t253 = -t236 * pkin(4) + pkin(6) * t291;
t277 = -t174 + t176 * t240 + pkin(1) + (qJ(3) - t253) * t238;
t276 = t174 - t205;
t186 = t295 * t242;
t219 = t242 * pkin(1) - qJ(2) * t244;
t275 = -t186 - t219;
t187 = t295 * t244;
t221 = pkin(1) * t244 + t242 * qJ(2);
t274 = -t187 - t221;
t273 = qJD(3) * t238;
t272 = V_base(4) * t219 + V_base(3);
t271 = V_base(5) * pkin(5) + V_base(1);
t183 = (-pkin(2) + t206) * t238;
t268 = -t183 + t292;
t136 = t242 * t276 + t244 * t278;
t267 = -t136 + t275;
t137 = -t242 * t278 + t244 * t276;
t266 = -t137 + t274;
t231 = V_base(6) + qJD(1);
t265 = t237 * t291;
t261 = qJD(2) * t242 + t271;
t259 = rSges(3,1) * t240 - rSges(3,2) * t238;
t258 = Icges(3,1) * t240 - t289;
t257 = -Icges(3,2) * t238 + t288;
t256 = Icges(3,5) * t240 - Icges(3,6) * t238;
t255 = -qJD(2) * t244 + t231 * t221 + V_base(2);
t254 = V_base(5) * t211 + t244 * t273 + t261;
t252 = -qJD(3) * t240 + V_base(4) * t186 + t272;
t251 = -t207 * t237 - t204 + t260;
t250 = t231 * t187 + t242 * t273 + t255;
t249 = qJD(4) * t196 + V_base(5) * t183 + t254;
t248 = qJD(4) * t287 + V_base(4) * t136 + t252;
t247 = qJD(4) * t194 + t231 * t137 + t250;
t246 = (-Icges(3,3) * t244 + t242 * t256) * V_base(5) + (Icges(3,3) * t242 + t244 * t256) * V_base(4) + (Icges(3,5) * t238 + Icges(3,6) * t240) * t231;
t179 = -Icges(3,6) * t244 + t242 * t257;
t180 = Icges(3,6) * t242 + t244 * t257;
t181 = -Icges(3,5) * t244 + t242 * t258;
t182 = Icges(3,5) * t242 + t244 * t258;
t209 = Icges(3,2) * t240 + t289;
t210 = Icges(3,1) * t238 + t288;
t245 = (-t180 * t238 + t182 * t240) * V_base(4) + (-t179 * t238 + t181 * t240) * V_base(5) + (-t209 * t238 + t210 * t240) * t231;
t234 = Icges(2,4) * t244;
t222 = rSges(2,1) * t244 - t242 * rSges(2,2);
t220 = t242 * rSges(2,1) + rSges(2,2) * t244;
t218 = Icges(2,1) * t244 - t290;
t217 = Icges(2,1) * t242 + t234;
t216 = -Icges(2,2) * t242 + t234;
t215 = Icges(2,2) * t244 + t290;
t214 = Icges(2,5) * t244 - Icges(2,6) * t242;
t213 = Icges(2,5) * t242 + Icges(2,6) * t244;
t212 = rSges(3,1) * t238 + rSges(3,2) * t240;
t203 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t202 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t201 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t193 = t243 * t239 + t241 * t265;
t192 = t239 * t263 + t284;
t185 = t242 * rSges(3,3) + t244 * t259;
t184 = -rSges(3,3) * t244 + t242 * t259;
t175 = -rSges(4,3) * t240 + (rSges(4,1) * t239 - rSges(4,2) * t237) * t238;
t173 = -Icges(4,5) * t240 + (Icges(4,1) * t239 - Icges(4,4) * t237) * t238;
t171 = -Icges(4,3) * t240 + (Icges(4,5) * t239 - Icges(4,6) * t237) * t238;
t168 = qJD(5) * t190 + t231;
t166 = V_base(5) * rSges(2,3) - t220 * t231 + t271;
t165 = t222 * t231 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t164 = t220 * V_base(4) - t222 * V_base(5) + V_base(3);
t159 = -t192 * t241 + t240 * t285;
t158 = -t191 * t241 + t238 * t285;
t157 = t191 * t243 + t238 * t286;
t156 = qJD(5) * t162 + V_base(4);
t155 = qJD(5) * t160 + V_base(5);
t153 = t197 * rSges(4,1) - t196 * rSges(4,2) + rSges(4,3) * t282;
t152 = rSges(4,1) * t195 - rSges(4,2) * t194 + rSges(4,3) * t283;
t151 = Icges(4,1) * t197 - Icges(4,4) * t196 + Icges(4,5) * t282;
t150 = Icges(4,1) * t195 - Icges(4,4) * t194 + Icges(4,5) * t283;
t147 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t282;
t146 = Icges(4,5) * t195 - Icges(4,6) * t194 + Icges(4,3) * t283;
t145 = rSges(5,1) * t191 - rSges(5,2) * t190 + rSges(5,3) * t287;
t144 = Icges(5,1) * t191 - Icges(5,4) * t190 + Icges(5,5) * t287;
t143 = Icges(5,4) * t191 - Icges(5,2) * t190 + Icges(5,6) * t287;
t141 = t159 * t244 - t242 * t193;
t140 = t159 * t242 + t193 * t244;
t139 = (t192 * t243 + t240 * t286) * t244 + t242 * (-t241 * t239 + t243 * t265);
t138 = (t192 * t242 - t237 * t262) * t243 + t194 * t241;
t135 = t253 * t240 + (-t206 + t176) * t238;
t132 = t212 * V_base(5) + (-t184 - t219) * t231 + t261;
t131 = t231 * t185 + (-pkin(5) - t212) * V_base(4) + t255;
t130 = t184 * V_base(4) + (-t185 - t221) * V_base(5) + t272;
t129 = rSges(5,1) * t163 - rSges(5,2) * t162 + rSges(5,3) * t196;
t128 = rSges(5,1) * t161 - rSges(5,2) * t160 + rSges(5,3) * t194;
t127 = Icges(5,1) * t163 - Icges(5,4) * t162 + Icges(5,5) * t196;
t126 = Icges(5,1) * t161 - Icges(5,4) * t160 + Icges(5,5) * t194;
t125 = Icges(5,4) * t163 - Icges(5,2) * t162 + Icges(5,6) * t196;
t124 = Icges(5,4) * t161 - Icges(5,2) * t160 + Icges(5,6) * t194;
t121 = rSges(6,1) * t157 + rSges(6,2) * t158 + rSges(6,3) * t190;
t120 = Icges(6,1) * t157 + Icges(6,4) * t158 + Icges(6,5) * t190;
t119 = Icges(6,4) * t157 + Icges(6,2) * t158 + Icges(6,6) * t190;
t118 = Icges(6,5) * t157 + Icges(6,6) * t158 + Icges(6,3) * t190;
t117 = -t242 * t251 + t244 * t277;
t116 = t242 * t277 + t244 * t251;
t115 = rSges(6,1) * t139 + rSges(6,2) * t141 + rSges(6,3) * t162;
t114 = rSges(6,1) * t138 + rSges(6,2) * t140 + rSges(6,3) * t160;
t113 = Icges(6,1) * t139 + Icges(6,4) * t141 + Icges(6,5) * t162;
t112 = Icges(6,1) * t138 + Icges(6,4) * t140 + Icges(6,5) * t160;
t111 = Icges(6,4) * t139 + Icges(6,2) * t141 + Icges(6,6) * t162;
t110 = Icges(6,4) * t138 + Icges(6,2) * t140 + Icges(6,6) * t160;
t109 = Icges(6,5) * t139 + Icges(6,6) * t141 + Icges(6,3) * t162;
t108 = Icges(6,5) * t138 + Icges(6,6) * t140 + Icges(6,3) * t160;
t107 = t175 * V_base(5) + (-t152 + t275) * t231 + t254;
t106 = t231 * t153 + (-t175 + t292) * V_base(4) + t250;
t105 = t152 * V_base(4) + (-t153 + t274) * V_base(5) + t252;
t104 = t145 * V_base(5) + (-t128 + t267) * t231 + t249;
t103 = t231 * t129 + (-t145 + t268) * V_base(4) + t247;
t102 = t128 * V_base(4) + (-t129 + t266) * V_base(5) + t248;
t101 = -t114 * t168 + t121 * t155 + t135 * V_base(5) + (-t116 + t267) * t231 + t249;
t100 = t168 * t115 + t231 * t117 - t156 * t121 + (-t135 + t268) * V_base(4) + t247;
t99 = t114 * t156 - t115 * t155 + t116 * V_base(4) + (-t117 + t266) * V_base(5) + t248;
t1 = m(1) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(2) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + t156 * ((t109 * t162 + t111 * t141 + t113 * t139) * t156 + (t108 * t162 + t110 * t141 + t112 * t139) * t155 + (t118 * t162 + t119 * t141 + t120 * t139) * t168) / 0.2e1 + t155 * ((t109 * t160 + t111 * t140 + t113 * t138) * t156 + (t108 * t160 + t110 * t140 + t112 * t138) * t155 + (t118 * t160 + t119 * t140 + t120 * t138) * t168) / 0.2e1 + t168 * ((t109 * t190 + t111 * t158 + t113 * t157) * t156 + (t108 * t190 + t110 * t158 + t112 * t157) * t155 + (t118 * t190 + t119 * t158 + t120 * t157) * t168) / 0.2e1 + ((t122 * t287 - t124 * t190 + t126 * t191 + t213 + (t179 - t146) * t240 + (-t148 * t237 + t150 * t239 + t181) * t238) * V_base(5) + (t123 * t287 - t125 * t190 + t127 * t191 + t214 + (t180 - t147) * t240 + (-t149 * t237 + t151 * t239 + t182) * t238) * V_base(4) + (t142 * t287 - t143 * t190 + t144 * t191 + Icges(2,3) + (t209 - t171) * t240 + (-t172 * t237 + t173 * t239 + t210) * t238) * t231) * t231 / 0.2e1 + (t242 * t246 + t244 * t245 + (-t143 * t162 + t144 * t163 + t171 * t282 + t197 * t173 + t196 * t297 + t214) * t231 + (-t124 * t162 + t126 * t163 + t146 * t282 + t197 * t150 + t196 * t299 - t242 * t215 + t217 * t244 + Icges(1,4)) * V_base(5) + (-t125 * t162 + t127 * t163 + t147 * t282 + t197 * t151 + t196 * t298 - t242 * t216 + t218 * t244 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t242 * t245 - t244 * t246 + (-t143 * t160 + t144 * t161 + t171 * t283 + t173 * t195 + t194 * t297 + t213) * t231 + (-t124 * t160 + t126 * t161 + t146 * t283 + t150 * t195 + t194 * t299 + t215 * t244 + t242 * t217 + Icges(1,2)) * V_base(5) + (-t125 * t160 + t127 * t161 + t147 * t283 + t151 * t195 + t194 * t298 + t216 * t244 + t242 * t218 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
