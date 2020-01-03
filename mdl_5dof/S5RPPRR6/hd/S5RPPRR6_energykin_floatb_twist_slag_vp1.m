% Calculate kinetic energy for
% S5RPPRR6
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:38
% EndTime: 2019-12-31 17:57:39
% DurationCPUTime: 1.63s
% Computational Cost: add. (1325->270), mult. (1130->377), div. (0->0), fcn. (966->10), ass. (0->138)
t180 = qJ(1) + pkin(8);
t172 = sin(t180);
t174 = cos(t180);
t185 = sin(qJ(1));
t187 = cos(qJ(1));
t237 = Icges(2,5) * t185 + Icges(3,5) * t172 + Icges(2,6) * t187 + Icges(3,6) * t174;
t236 = Icges(2,5) * t187 + Icges(3,5) * t174 - Icges(2,6) * t185 - Icges(3,6) * t172;
t234 = pkin(1) * t185;
t233 = pkin(1) * t187;
t181 = sin(pkin(9));
t232 = pkin(3) * t181;
t182 = cos(pkin(9));
t231 = pkin(3) * t182;
t230 = -pkin(5) - qJ(2);
t229 = Icges(2,4) * t185;
t228 = Icges(3,4) * t172;
t227 = Icges(4,4) * t181;
t226 = Icges(4,4) * t182;
t179 = pkin(9) + qJ(4);
t171 = sin(t179);
t225 = Icges(5,4) * t171;
t173 = cos(t179);
t224 = Icges(5,4) * t173;
t223 = t171 * t172;
t222 = t171 * t174;
t184 = sin(qJ(5));
t221 = t172 * t184;
t186 = cos(qJ(5));
t220 = t172 * t186;
t219 = t174 * t184;
t218 = t174 * t186;
t216 = qJD(5) * t171;
t175 = V_base(6) + qJD(1);
t215 = t175 * t233 + V_base(2);
t214 = V_base(5) * pkin(5) + V_base(1);
t154 = qJD(4) * t172 + V_base(4);
t141 = pkin(2) * t172 - qJ(3) * t174;
t211 = -t141 - t234;
t143 = pkin(2) * t174 + qJ(3) * t172;
t210 = -t143 - t233;
t209 = V_base(5) * qJ(2) + t214;
t208 = V_base(4) * t234 + qJD(2) + V_base(3);
t96 = -pkin(6) * t174 + t172 * t231;
t207 = t211 - t96;
t206 = qJD(3) * t172 + t209;
t205 = pkin(4) * t173 + pkin(7) * t171;
t204 = V_base(4) * t141 + t208;
t153 = -qJD(4) * t174 + V_base(5);
t203 = rSges(4,1) * t182 - rSges(4,2) * t181;
t202 = rSges(5,1) * t173 - rSges(5,2) * t171;
t201 = Icges(4,1) * t182 - t227;
t200 = Icges(5,1) * t173 - t225;
t199 = -Icges(4,2) * t181 + t226;
t198 = -Icges(5,2) * t171 + t224;
t197 = Icges(4,5) * t182 - Icges(4,6) * t181;
t196 = Icges(5,5) * t173 - Icges(5,6) * t171;
t195 = V_base(5) * t232 + t206;
t194 = -qJD(3) * t174 + t175 * t143 + t215;
t193 = (Icges(5,5) * t171 + Icges(5,6) * t173) * t175 + (-Icges(5,3) * t174 + t172 * t196) * t153 + (Icges(5,3) * t172 + t174 * t196) * t154;
t192 = (-Icges(4,3) * t174 + t172 * t197) * V_base(5) + (Icges(4,3) * t172 + t174 * t197) * V_base(4) + (Icges(4,5) * t181 + Icges(4,6) * t182) * t175;
t97 = pkin(6) * t172 + t174 * t231;
t191 = V_base(4) * t96 + (t210 - t97) * V_base(5) + t204;
t190 = t175 * t97 + (t230 - t232) * V_base(4) + t194;
t100 = -Icges(5,6) * t174 + t172 * t198;
t101 = Icges(5,6) * t172 + t174 * t198;
t102 = -Icges(5,5) * t174 + t172 * t200;
t103 = Icges(5,5) * t172 + t174 * t200;
t134 = Icges(5,2) * t173 + t225;
t137 = Icges(5,1) * t171 + t224;
t189 = (-t101 * t171 + t103 * t173) * t154 + (-t100 * t171 + t102 * t173) * t153 + (-t134 * t171 + t137 * t173) * t175;
t109 = -Icges(4,6) * t174 + t172 * t199;
t110 = Icges(4,6) * t172 + t174 * t199;
t111 = -Icges(4,5) * t174 + t172 * t201;
t112 = Icges(4,5) * t172 + t174 * t201;
t151 = Icges(4,2) * t182 + t227;
t152 = Icges(4,1) * t181 + t226;
t188 = (-t110 * t181 + t112 * t182) * V_base(4) + (-t109 * t181 + t111 * t182) * V_base(5) + (-t151 * t181 + t152 * t182) * t175;
t177 = Icges(2,4) * t187;
t169 = Icges(3,4) * t174;
t163 = rSges(2,1) * t187 - t185 * rSges(2,2);
t162 = t185 * rSges(2,1) + rSges(2,2) * t187;
t161 = Icges(2,1) * t187 - t229;
t160 = Icges(2,1) * t185 + t177;
t159 = -Icges(2,2) * t185 + t177;
t158 = Icges(2,2) * t187 + t229;
t155 = rSges(4,1) * t181 + rSges(4,2) * t182;
t149 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t148 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t147 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t146 = -qJD(5) * t173 + t175;
t145 = pkin(4) * t171 - pkin(7) * t173;
t144 = rSges(3,1) * t174 - rSges(3,2) * t172;
t142 = rSges(3,1) * t172 + rSges(3,2) * t174;
t140 = rSges(5,1) * t171 + rSges(5,2) * t173;
t139 = Icges(3,1) * t174 - t228;
t138 = Icges(3,1) * t172 + t169;
t136 = -Icges(3,2) * t172 + t169;
t135 = Icges(3,2) * t174 + t228;
t128 = t173 * t218 + t221;
t127 = -t173 * t219 + t220;
t126 = t173 * t220 - t219;
t125 = -t173 * t221 - t218;
t124 = t174 * t216 + t154;
t123 = t172 * t216 + t153;
t122 = t205 * t174;
t121 = t205 * t172;
t120 = -rSges(6,3) * t173 + (rSges(6,1) * t186 - rSges(6,2) * t184) * t171;
t119 = V_base(5) * rSges(2,3) - t162 * t175 + t214;
t118 = t163 * t175 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t117 = -Icges(6,5) * t173 + (Icges(6,1) * t186 - Icges(6,4) * t184) * t171;
t116 = -Icges(6,6) * t173 + (Icges(6,4) * t186 - Icges(6,2) * t184) * t171;
t115 = -Icges(6,3) * t173 + (Icges(6,5) * t186 - Icges(6,6) * t184) * t171;
t114 = rSges(4,3) * t172 + t174 * t203;
t113 = -rSges(4,3) * t174 + t172 * t203;
t106 = t162 * V_base(4) - t163 * V_base(5) + V_base(3);
t105 = rSges(5,3) * t172 + t174 * t202;
t104 = -rSges(5,3) * t174 + t172 * t202;
t93 = V_base(5) * rSges(3,3) + (-t142 - t234) * t175 + t209;
t92 = t144 * t175 + (-rSges(3,3) + t230) * V_base(4) + t215;
t91 = V_base(4) * t142 + (-t144 - t233) * V_base(5) + t208;
t90 = rSges(6,1) * t128 + rSges(6,2) * t127 + rSges(6,3) * t222;
t89 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t223;
t88 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t222;
t87 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t223;
t86 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t222;
t85 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t223;
t84 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t222;
t83 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t223;
t82 = t155 * V_base(5) + (-t113 + t211) * t175 + t206;
t81 = t114 * t175 + (-t155 + t230) * V_base(4) + t194;
t80 = V_base(4) * t113 + (-t114 + t210) * V_base(5) + t204;
t79 = t140 * t153 + (-t104 + t207) * t175 + t195;
t78 = t105 * t175 - t140 * t154 + t190;
t77 = t154 * t104 - t153 * t105 + t191;
t76 = t120 * t123 + t145 * t153 - t146 * t89 + (-t121 + t207) * t175 + t195;
t75 = -t120 * t124 + t122 * t175 - t145 * t154 + t146 * t90 + t190;
t74 = t154 * t121 - t153 * t122 - t123 * t90 + t124 * t89 + t191;
t1 = m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(2) * (t106 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + t154 * (t172 * t193 + t174 * t189) / 0.2e1 + t153 * (t172 * t189 - t174 * t193) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t124 * ((t127 * t86 + t128 * t88 + t84 * t222) * t124 + (t127 * t85 + t128 * t87 + t222 * t83) * t123 + (t115 * t222 + t116 * t127 + t117 * t128) * t146) / 0.2e1 + t123 * ((t125 * t86 + t126 * t88 + t223 * t84) * t124 + (t125 * t85 + t126 * t87 + t83 * t223) * t123 + (t115 * t223 + t116 * t125 + t117 * t126) * t146) / 0.2e1 + t146 * ((-t115 * t146 - t83 * t123 - t84 * t124) * t173 + ((-t184 * t86 + t186 * t88) * t124 + (-t184 * t85 + t186 * t87) * t123 + (-t116 * t184 + t117 * t186) * t146) * t171) / 0.2e1 + ((t101 * t173 + t103 * t171) * t154 + (t100 * t173 + t102 * t171) * t153 + (t109 * t182 + t111 * t181 + t237) * V_base(5) + (t110 * t182 + t112 * t181 + t236) * V_base(4) + (t134 * t173 + t137 * t171 + t151 * t182 + t152 * t181 + Icges(2,3) + Icges(3,3)) * t175) * t175 / 0.2e1 + (t172 * t192 + t174 * t188 + t236 * t175 + (-t135 * t172 + t138 * t174 - t185 * t158 + t160 * t187 + Icges(1,4)) * V_base(5) + (-t136 * t172 + t139 * t174 - t185 * t159 + t161 * t187 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t172 * t188 - t174 * t192 + t237 * t175 + (t135 * t174 + t138 * t172 + t158 * t187 + t185 * t160 + Icges(1,2)) * V_base(5) + (t136 * t174 + t139 * t172 + t159 * t187 + t185 * t161 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
