% Calculate kinetic energy for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:21
% EndTime: 2020-01-03 11:30:24
% DurationCPUTime: 3.26s
% Computational Cost: add. (1399->337), mult. (1865->476), div. (0->0), fcn. (1815->10), ass. (0->156)
t221 = sin(qJ(1));
t222 = cos(qJ(1));
t217 = sin(pkin(8));
t219 = cos(pkin(8));
t263 = Icges(3,4) * t219;
t236 = -Icges(3,2) * t217 + t263;
t156 = -Icges(3,6) * t222 + t221 * t236;
t264 = Icges(3,4) * t217;
t237 = Icges(3,1) * t219 - t264;
t158 = -Icges(3,5) * t222 + t221 * t237;
t265 = Icges(2,4) * t222;
t273 = Icges(2,1) * t221 - t156 * t217 + t158 * t219 + t265;
t157 = -Icges(3,6) * t221 - t222 * t236;
t159 = -Icges(3,5) * t221 - t222 * t237;
t212 = Icges(2,4) * t221;
t272 = -Icges(2,1) * t222 - t157 * t217 + t159 * t219 + t212;
t215 = pkin(9) + qJ(4);
t208 = cos(t215);
t249 = pkin(4) * t208;
t271 = pkin(7) * t217 + t219 * t249;
t218 = cos(pkin(9));
t267 = t218 * pkin(3);
t270 = pkin(6) * t217 + t219 * t267;
t189 = pkin(2) * t217 - qJ(3) * t219;
t266 = -pkin(5) - t189;
t216 = sin(pkin(9));
t262 = t216 * t222;
t261 = t217 * t221;
t260 = t217 * t222;
t259 = t219 * t222;
t209 = qJ(5) + t215;
t204 = sin(t209);
t258 = t221 * t204;
t205 = cos(t209);
t257 = t221 * t205;
t207 = sin(t215);
t256 = t221 * t207;
t255 = t221 * t208;
t254 = t221 * t216;
t253 = t221 * t218;
t238 = pkin(2) * t219 + qJ(3) * t217;
t173 = t238 * t221;
t198 = t221 * pkin(1) - qJ(2) * t222;
t251 = -t173 - t198;
t174 = t238 * t222;
t200 = -pkin(1) * t222 - t221 * qJ(2);
t250 = t174 - t200;
t247 = qJD(3) * t217;
t246 = qJD(4) * t217;
t245 = -qJD(4) - qJD(5);
t244 = V_base(5) * t200 + V_base(1);
t243 = V_base(6) * pkin(5) + V_base(2);
t185 = t221 * t246 + V_base(5);
t210 = V_base(4) + qJD(1);
t240 = pkin(4) * t207;
t239 = rSges(3,1) * t219 - rSges(3,2) * t217;
t235 = Icges(3,5) * t219 - Icges(3,6) * t217;
t187 = Icges(3,2) * t219 + t264;
t188 = Icges(3,1) * t217 + t263;
t232 = t187 * t217 - t188 * t219;
t231 = -qJD(2) * t221 + t210 * t198 + V_base(3);
t230 = -qJD(2) * t222 + t243;
t229 = -qJD(3) * t219 - V_base(5) * t174 + t244;
t228 = V_base(6) * t189 + t221 * t247 + t230;
t227 = -(-Icges(3,3) * t222 + t221 * t235) * V_base(5) - (-Icges(3,3) * t221 - t222 * t235) * V_base(6) - (Icges(3,5) * t217 + Icges(3,6) * t219) * t210;
t226 = t210 * t173 - t222 * t247 + t231;
t128 = -pkin(3) * t254 - t222 * t270;
t138 = -pkin(6) * t219 + t217 * t267;
t225 = V_base(6) * t138 + (-t128 + t250) * t210 + t228;
t127 = -pkin(3) * t262 + t221 * t270;
t224 = V_base(5) * t128 + (-t127 + t251) * V_base(6) + t229;
t223 = t210 * t127 + (-t138 + t266) * V_base(5) + t226;
t201 = -rSges(2,1) * t222 + t221 * rSges(2,2);
t199 = t221 * rSges(2,1) + rSges(2,2) * t222;
t195 = Icges(2,2) * t221 - t265;
t194 = Icges(2,2) * t222 + t212;
t193 = -Icges(2,5) * t222 + Icges(2,6) * t221;
t192 = Icges(2,5) * t221 + Icges(2,6) * t222;
t191 = -qJD(4) * t219 + t210;
t190 = rSges(3,1) * t217 + rSges(3,2) * t219;
t184 = -t222 * t246 + V_base(6);
t183 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t182 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t181 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t175 = t219 * t245 + t210;
t172 = -t218 * t259 - t254;
t171 = t216 * t259 - t253;
t170 = t219 * t253 - t262;
t169 = -t218 * t222 - t219 * t254;
t167 = qJD(5) * t261 + t185;
t166 = t245 * t260 + V_base(6);
t165 = -t221 * rSges(3,3) - t222 * t239;
t164 = -rSges(3,3) * t222 + t221 * t239;
t163 = -t208 * t259 - t256;
t162 = t207 * t259 - t255;
t161 = -t207 * t222 + t219 * t255;
t160 = -t208 * t222 - t219 * t256;
t152 = -rSges(4,3) * t219 + (rSges(4,1) * t218 - rSges(4,2) * t216) * t217;
t151 = -Icges(4,5) * t219 + (Icges(4,1) * t218 - Icges(4,4) * t216) * t217;
t150 = -Icges(4,6) * t219 + (Icges(4,4) * t218 - Icges(4,2) * t216) * t217;
t149 = -Icges(4,3) * t219 + (Icges(4,5) * t218 - Icges(4,6) * t216) * t217;
t148 = -t205 * t259 - t258;
t147 = t204 * t259 - t257;
t146 = -t204 * t222 + t219 * t257;
t145 = -t205 * t222 - t219 * t258;
t144 = -rSges(5,3) * t219 + (rSges(5,1) * t208 - rSges(5,2) * t207) * t217;
t143 = -Icges(5,5) * t219 + (Icges(5,1) * t208 - Icges(5,4) * t207) * t217;
t142 = -Icges(5,6) * t219 + (Icges(5,4) * t208 - Icges(5,2) * t207) * t217;
t141 = -Icges(5,3) * t219 + (Icges(5,5) * t208 - Icges(5,6) * t207) * t217;
t140 = V_base(6) * rSges(2,3) - t201 * t210 + t243;
t139 = t199 * t210 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t137 = -rSges(6,3) * t219 + (rSges(6,1) * t205 - rSges(6,2) * t204) * t217;
t136 = -Icges(6,5) * t219 + (Icges(6,1) * t205 - Icges(6,4) * t204) * t217;
t135 = -Icges(6,6) * t219 + (Icges(6,4) * t205 - Icges(6,2) * t204) * t217;
t134 = -Icges(6,3) * t219 + (Icges(6,5) * t205 - Icges(6,6) * t204) * t217;
t133 = -t199 * V_base(6) + t201 * V_base(5) + V_base(1);
t131 = -pkin(7) * t219 + t217 * t249;
t130 = t172 * rSges(4,1) + t171 * rSges(4,2) - rSges(4,3) * t260;
t129 = rSges(4,1) * t170 + rSges(4,2) * t169 + rSges(4,3) * t261;
t126 = Icges(4,1) * t172 + Icges(4,4) * t171 - Icges(4,5) * t260;
t125 = Icges(4,1) * t170 + Icges(4,4) * t169 + Icges(4,5) * t261;
t124 = Icges(4,4) * t172 + Icges(4,2) * t171 - Icges(4,6) * t260;
t123 = Icges(4,4) * t170 + Icges(4,2) * t169 + Icges(4,6) * t261;
t122 = Icges(4,5) * t172 + Icges(4,6) * t171 - Icges(4,3) * t260;
t121 = Icges(4,5) * t170 + Icges(4,6) * t169 + Icges(4,3) * t261;
t118 = t163 * rSges(5,1) + t162 * rSges(5,2) - rSges(5,3) * t260;
t117 = rSges(5,1) * t161 + rSges(5,2) * t160 + rSges(5,3) * t261;
t116 = Icges(5,1) * t163 + Icges(5,4) * t162 - Icges(5,5) * t260;
t115 = Icges(5,1) * t161 + Icges(5,4) * t160 + Icges(5,5) * t261;
t114 = Icges(5,4) * t163 + Icges(5,2) * t162 - Icges(5,6) * t260;
t113 = Icges(5,4) * t161 + Icges(5,2) * t160 + Icges(5,6) * t261;
t112 = Icges(5,5) * t163 + Icges(5,6) * t162 - Icges(5,3) * t260;
t111 = Icges(5,5) * t161 + Icges(5,6) * t160 + Icges(5,3) * t261;
t110 = t148 * rSges(6,1) + t147 * rSges(6,2) - rSges(6,3) * t260;
t109 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t261;
t108 = Icges(6,1) * t148 + Icges(6,4) * t147 - Icges(6,5) * t260;
t107 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t261;
t106 = Icges(6,4) * t148 + Icges(6,2) * t147 - Icges(6,6) * t260;
t105 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t261;
t104 = Icges(6,5) * t148 + Icges(6,6) * t147 - Icges(6,3) * t260;
t103 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t261;
t102 = V_base(6) * t190 + (-t165 - t200) * t210 + t230;
t101 = t164 * t210 + (-pkin(5) - t190) * V_base(5) + t231;
t100 = -t240 * t221 - t222 * t271;
t99 = t221 * t271 - t240 * t222;
t98 = t165 * V_base(5) + (-t164 - t198) * V_base(6) + t244;
t97 = V_base(6) * t152 + (-t130 + t250) * t210 + t228;
t96 = t210 * t129 + (-t152 + t266) * V_base(5) + t226;
t95 = t130 * V_base(5) + (-t129 + t251) * V_base(6) + t229;
t94 = -t191 * t118 + t184 * t144 + t225;
t93 = t191 * t117 - t185 * t144 + t223;
t92 = -t117 * t184 + t118 * t185 + t224;
t91 = -t191 * t100 - t175 * t110 + t184 * t131 + t166 * t137 + t225;
t90 = t175 * t109 - t185 * t131 - t167 * t137 + t191 * t99 + t223;
t89 = t100 * t185 - t109 * t166 + t110 * t167 - t184 * t99 + t224;
t1 = m(1) * (t181 ^ 2 + t182 ^ 2 + t183 ^ 2) / 0.2e1 + m(2) * (t133 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(3) * (t101 ^ 2 + t102 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t191 * ((-t111 * t185 - t112 * t184 - t141 * t191) * t219 + ((-t142 * t207 + t143 * t208) * t191 + (-t113 * t207 + t115 * t208) * t185 + (-t114 * t207 + t116 * t208) * t184) * t217) / 0.2e1 + t185 * ((t141 * t261 + t142 * t160 + t143 * t161) * t191 + (t111 * t261 + t160 * t113 + t161 * t115) * t185 + (t112 * t261 + t114 * t160 + t116 * t161) * t184) / 0.2e1 + t184 * ((-t141 * t260 + t162 * t142 + t163 * t143) * t191 + (-t111 * t260 + t162 * t113 + t163 * t115) * t185 + (-t112 * t260 + t162 * t114 + t163 * t116) * t184) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t175 * ((-t103 * t167 - t104 * t166 - t134 * t175) * t219 + ((-t135 * t204 + t136 * t205) * t175 + (-t105 * t204 + t107 * t205) * t167 + (-t106 * t204 + t108 * t205) * t166) * t217) / 0.2e1 + t167 * ((t134 * t261 + t135 * t145 + t136 * t146) * t175 + (t103 * t261 + t145 * t105 + t146 * t107) * t167 + (t104 * t261 + t106 * t145 + t108 * t146) * t166) / 0.2e1 + t166 * ((-t134 * t260 + t147 * t135 + t148 * t136) * t175 + (-t103 * t260 + t147 * t105 + t148 * t107) * t167 + (-t104 * t260 + t147 * t106 + t148 * t108) * t166) / 0.2e1 + ((t193 + (t157 - t122) * t219 + (-t124 * t216 + t126 * t218 + t159) * t217) * V_base(6) + (t192 + (t156 - t121) * t219 + (-t123 * t216 + t125 * t218 + t158) * t217) * V_base(5) + (Icges(2,3) + (t187 - t149) * t219 + (-t150 * t216 + t151 * t218 + t188) * t217) * t210) * t210 / 0.2e1 + (t227 * t222 + (t149 * t261 + t150 * t169 + t151 * t170 - t232 * t221 + t192) * t210 + (t122 * t261 + t124 * t169 + t126 * t170 + t195 * t222 + t221 * t272 + Icges(1,6)) * V_base(6) + (t121 * t261 + t169 * t123 + t170 * t125 + t222 * t194 + t273 * t221 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t227 * t221 + (-t149 * t260 + t171 * t150 + t172 * t151 + t232 * t222 + t193) * t210 + (-t122 * t260 + t171 * t124 + t172 * t126 + t221 * t195 - t272 * t222 + Icges(1,3)) * V_base(6) + (-t121 * t260 + t171 * t123 + t172 * t125 + t221 * t194 - t222 * t273 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
