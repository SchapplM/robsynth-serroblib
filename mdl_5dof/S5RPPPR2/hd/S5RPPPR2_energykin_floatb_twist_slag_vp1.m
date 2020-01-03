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
% m_mdh [6x1]
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:16
% EndTime: 2020-01-03 11:22:19
% DurationCPUTime: 2.68s
% Computational Cost: add. (1293->336), mult. (2923->444), div. (0->0), fcn. (3349->10), ass. (0->151)
t231 = sin(pkin(8));
t234 = cos(pkin(7));
t238 = cos(qJ(1));
t233 = cos(pkin(8));
t236 = sin(qJ(1));
t267 = t236 * t233;
t196 = -t231 * t238 + t234 * t267;
t230 = sin(pkin(9));
t232 = sin(pkin(7));
t276 = cos(pkin(9));
t256 = t232 * t276;
t168 = t196 * t230 - t236 * t256;
t271 = t232 * t236;
t169 = t196 * t276 + t230 * t271;
t268 = t236 * t231;
t195 = t233 * t238 + t234 * t268;
t130 = Icges(5,5) * t169 - Icges(5,6) * t168 + Icges(5,3) * t195;
t153 = Icges(4,4) * t196 - Icges(4,2) * t195 + Icges(4,6) * t271;
t283 = t130 - t153;
t269 = t234 * t238;
t198 = -t233 * t269 - t268;
t170 = t198 * t230 + t238 * t256;
t270 = t232 * t238;
t171 = t198 * t276 - t230 * t270;
t197 = t231 * t269 - t267;
t131 = Icges(5,5) * t171 - Icges(5,6) * t170 - Icges(5,3) * t197;
t154 = Icges(4,4) * t198 + Icges(4,2) * t197 - Icges(4,6) * t270;
t282 = t131 - t154;
t193 = t232 * t233 * t230 + t234 * t276;
t194 = -t234 * t230 + t233 * t256;
t272 = t231 * t232;
t147 = Icges(5,5) * t194 - Icges(5,6) * t193 + Icges(5,3) * t272;
t177 = -Icges(4,6) * t234 + (Icges(4,4) * t233 - Icges(4,2) * t231) * t232;
t281 = t147 - t177;
t273 = Icges(3,4) * t234;
t252 = -Icges(3,2) * t232 + t273;
t183 = -Icges(3,6) * t238 + t236 * t252;
t274 = Icges(3,4) * t232;
t253 = Icges(3,1) * t234 - t274;
t185 = -Icges(3,5) * t238 + t236 * t253;
t275 = Icges(2,4) * t238;
t280 = Icges(2,1) * t236 - t183 * t232 + t185 * t234 + t275;
t184 = -Icges(3,6) * t236 - t238 * t252;
t186 = -Icges(3,5) * t236 - t238 * t253;
t228 = Icges(2,4) * t236;
t279 = -Icges(2,1) * t238 - t184 * t232 + t186 * t234 + t228;
t211 = pkin(2) * t232 - qJ(3) * t234;
t277 = -pkin(5) - t211;
t254 = pkin(2) * t234 + qJ(3) * t232;
t200 = t254 * t236;
t219 = t236 * pkin(1) - qJ(2) * t238;
t266 = -t200 - t219;
t201 = t254 * t238;
t221 = -pkin(1) * t238 - t236 * qJ(2);
t265 = t201 - t221;
t264 = qJD(3) * t232;
t263 = V_base(5) * t221 + V_base(1);
t262 = V_base(6) * pkin(5) + V_base(2);
t199 = (pkin(3) * t233 + qJ(4) * t231) * t232;
t259 = -t199 + t277;
t162 = pkin(3) * t196 + qJ(4) * t195;
t258 = -t162 + t266;
t163 = pkin(3) * t198 - qJ(4) * t197;
t257 = -t163 + t265;
t227 = V_base(4) + qJD(1);
t255 = rSges(3,1) * t234 - rSges(3,2) * t232;
t251 = Icges(3,5) * t234 - Icges(3,6) * t232;
t209 = Icges(3,2) * t234 + t274;
t210 = Icges(3,1) * t232 + t273;
t248 = t209 * t232 - t210 * t234;
t247 = -qJD(2) * t236 + t227 * t219 + V_base(3);
t246 = -qJD(2) * t238 + t262;
t245 = -qJD(3) * t234 - V_base(5) * t201 + t263;
t244 = V_base(6) * t211 + t236 * t264 + t246;
t243 = qJD(4) * t272 + V_base(5) * t163 + t245;
t242 = qJD(4) * t195 + V_base(6) * t199 + t244;
t241 = -(-Icges(3,3) * t238 + t236 * t251) * V_base(5) - (-Icges(3,3) * t236 - t238 * t251) * V_base(6) - (Icges(3,5) * t232 + Icges(3,6) * t234) * t227;
t240 = t227 * t200 - t238 * t264 + t247;
t239 = -qJD(4) * t197 + t227 * t162 + t240;
t237 = cos(qJ(5));
t235 = sin(qJ(5));
t222 = -rSges(2,1) * t238 + t236 * rSges(2,2);
t220 = t236 * rSges(2,1) + rSges(2,2) * t238;
t216 = Icges(2,2) * t236 - t275;
t215 = Icges(2,2) * t238 + t228;
t214 = -Icges(2,5) * t238 + Icges(2,6) * t236;
t213 = Icges(2,5) * t236 + Icges(2,6) * t238;
t212 = rSges(3,1) * t232 + rSges(3,2) * t234;
t207 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t206 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t205 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t188 = -t236 * rSges(3,3) - t238 * t255;
t187 = -rSges(3,3) * t238 + t236 * t255;
t179 = -rSges(4,3) * t234 + (rSges(4,1) * t233 - rSges(4,2) * t231) * t232;
t178 = -Icges(4,5) * t234 + (Icges(4,1) * t233 - Icges(4,4) * t231) * t232;
t176 = -Icges(4,3) * t234 + (Icges(4,5) * t233 - Icges(4,6) * t231) * t232;
t175 = qJD(5) * t193 + t227;
t174 = V_base(6) * rSges(2,3) - t222 * t227 + t262;
t173 = t220 * t227 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t172 = -t220 * V_base(6) + t222 * V_base(5) + V_base(1);
t167 = t194 * t237 + t235 * t272;
t166 = -t194 * t235 + t237 * t272;
t165 = qJD(5) * t168 + V_base(5);
t164 = qJD(5) * t170 + V_base(6);
t161 = pkin(4) * t194 + pkin(6) * t193;
t158 = t198 * rSges(4,1) + t197 * rSges(4,2) - rSges(4,3) * t270;
t157 = rSges(4,1) * t196 - rSges(4,2) * t195 + rSges(4,3) * t271;
t156 = Icges(4,1) * t198 + Icges(4,4) * t197 - Icges(4,5) * t270;
t155 = Icges(4,1) * t196 - Icges(4,4) * t195 + Icges(4,5) * t271;
t152 = Icges(4,5) * t198 + Icges(4,6) * t197 - Icges(4,3) * t270;
t151 = Icges(4,5) * t196 - Icges(4,6) * t195 + Icges(4,3) * t271;
t150 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t272;
t149 = Icges(5,1) * t194 - Icges(5,4) * t193 + Icges(5,5) * t272;
t148 = Icges(5,4) * t194 - Icges(5,2) * t193 + Icges(5,6) * t272;
t146 = t171 * t237 - t197 * t235;
t145 = -t171 * t235 - t197 * t237;
t144 = t169 * t237 + t195 * t235;
t143 = -t169 * t235 + t195 * t237;
t142 = pkin(4) * t171 + pkin(6) * t170;
t141 = pkin(4) * t169 + pkin(6) * t168;
t140 = V_base(6) * t212 + (-t188 - t221) * t227 + t246;
t139 = t187 * t227 + (-pkin(5) - t212) * V_base(5) + t247;
t138 = t188 * V_base(5) + (-t187 - t219) * V_base(6) + t263;
t137 = rSges(5,1) * t171 - rSges(5,2) * t170 - rSges(5,3) * t197;
t136 = rSges(5,1) * t169 - rSges(5,2) * t168 + rSges(5,3) * t195;
t135 = Icges(5,1) * t171 - Icges(5,4) * t170 - Icges(5,5) * t197;
t134 = Icges(5,1) * t169 - Icges(5,4) * t168 + Icges(5,5) * t195;
t133 = Icges(5,4) * t171 - Icges(5,2) * t170 - Icges(5,6) * t197;
t132 = Icges(5,4) * t169 - Icges(5,2) * t168 + Icges(5,6) * t195;
t129 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t193;
t128 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t193;
t127 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t193;
t126 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t193;
t125 = V_base(6) * t179 + (-t158 + t265) * t227 + t244;
t124 = t227 * t157 + (-t179 + t277) * V_base(5) + t240;
t123 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t170;
t122 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t168;
t121 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t170;
t120 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t168;
t119 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t170;
t118 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t168;
t117 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t170;
t116 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t168;
t115 = t158 * V_base(5) + (-t157 + t266) * V_base(6) + t245;
t114 = V_base(6) * t150 + (-t137 + t257) * t227 + t242;
t113 = t227 * t136 + (-t150 + t259) * V_base(5) + t239;
t112 = t137 * V_base(5) + (-t136 + t258) * V_base(6) + t243;
t111 = -t175 * t123 + t164 * t129 + V_base(6) * t161 + (-t142 + t257) * t227 + t242;
t110 = t175 * t122 - t165 * t129 + t227 * t141 + (-t161 + t259) * V_base(5) + t239;
t109 = -t122 * t164 + t123 * t165 + t142 * V_base(5) + (-t141 + t258) * V_base(6) + t243;
t1 = m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(2) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(3) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t175 * ((t126 * t193 + t127 * t166 + t128 * t167) * t175 + (t116 * t193 + t118 * t166 + t120 * t167) * t165 + (t117 * t193 + t119 * t166 + t121 * t167) * t164) / 0.2e1 + t165 * ((t126 * t168 + t127 * t143 + t128 * t144) * t175 + (t116 * t168 + t118 * t143 + t120 * t144) * t165 + (t117 * t168 + t119 * t143 + t121 * t144) * t164) / 0.2e1 + t164 * ((t126 * t170 + t127 * t145 + t128 * t146) * t175 + (t116 * t170 + t118 * t145 + t120 * t146) * t165 + (t117 * t170 + t119 * t145 + t121 * t146) * t164) / 0.2e1 + ((t131 * t272 - t133 * t193 + t135 * t194 + t214 + (t184 - t152) * t234 + (-t154 * t231 + t156 * t233 + t186) * t232) * V_base(6) + (t130 * t272 - t132 * t193 + t134 * t194 + t213 + (t183 - t151) * t234 + (-t153 * t231 + t155 * t233 + t185) * t232) * V_base(5) + (t147 * t272 - t148 * t193 + t149 * t194 + Icges(2,3) + (t209 - t176) * t234 + (-t177 * t231 + t178 * t233 + t210) * t232) * t227) * t227 / 0.2e1 + (t241 * t238 + (-t148 * t168 + t149 * t169 + t176 * t271 + t178 * t196 + t195 * t281 - t248 * t236 + t213) * t227 + (-t133 * t168 + t135 * t169 + t152 * t271 + t156 * t196 + t195 * t282 + t216 * t238 + t236 * t279 + Icges(1,6)) * V_base(6) + (-t132 * t168 + t134 * t169 + t151 * t271 + t155 * t196 + t283 * t195 + t215 * t238 + t280 * t236 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t241 * t236 + (-t148 * t170 + t149 * t171 - t176 * t270 + t198 * t178 - t197 * t281 + t248 * t238 + t214) * t227 + (-t133 * t170 + t135 * t171 - t152 * t270 + t198 * t156 - t282 * t197 + t236 * t216 - t279 * t238 + Icges(1,3)) * V_base(6) + (-t132 * t170 + t134 * t171 - t151 * t270 + t198 * t155 - t197 * t283 + t236 * t215 - t280 * t238 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
