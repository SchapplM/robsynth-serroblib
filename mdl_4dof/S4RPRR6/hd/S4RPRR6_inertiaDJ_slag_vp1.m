% Calculate time derivative of joint inertia matrix for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:34
% DurationCPUTime: 2.64s
% Computational Cost: add. (4275->287), mult. (4240->424), div. (0->0), fcn. (3290->8), ass. (0->170)
t120 = sin(qJ(1));
t121 = cos(qJ(1));
t113 = pkin(7) + qJ(3);
t106 = qJ(4) + t113;
t100 = cos(t106);
t204 = rSges(5,1) * t100;
t99 = sin(t106);
t206 = rSges(5,2) * t99;
t154 = t204 - t206;
t67 = t120 * rSges(5,3) + t121 * t154;
t114 = qJD(3) + qJD(4);
t202 = rSges(5,2) * t100;
t87 = rSges(5,1) * t99 + t202;
t132 = t87 * t114;
t104 = sin(t113);
t169 = qJD(3) * t104;
t164 = pkin(3) * t169;
t220 = t164 + t132;
t141 = Icges(5,5) * t100 - Icges(5,6) * t99;
t60 = -Icges(5,3) * t121 + t120 * t141;
t219 = qJD(1) * t60;
t105 = cos(t113);
t186 = Icges(4,4) * t105;
t138 = -Icges(4,2) * t104 + t186;
t71 = Icges(4,6) * t120 + t121 * t138;
t187 = Icges(4,4) * t104;
t140 = Icges(4,1) * t105 - t187;
t73 = Icges(4,5) * t120 + t121 * t140;
t146 = t104 * t71 - t105 * t73;
t218 = t120 * t146;
t185 = Icges(5,4) * t100;
t142 = -Icges(5,2) * t99 + t185;
t63 = Icges(5,6) * t120 + t121 * t142;
t201 = Icges(5,4) * t99;
t143 = Icges(5,1) * t100 - t201;
t65 = Icges(5,5) * t120 + t121 * t143;
t152 = t100 * t65 - t63 * t99;
t217 = t120 * t152;
t70 = -Icges(4,6) * t121 + t120 * t138;
t72 = -Icges(4,5) * t121 + t120 * t140;
t147 = t104 * t70 - t105 * t72;
t216 = t121 * t147;
t62 = -Icges(5,6) * t121 + t120 * t142;
t64 = -Icges(5,5) * t121 + t120 * t143;
t153 = t100 * t64 - t62 * t99;
t215 = t121 * t153;
t119 = -pkin(5) - qJ(2);
t112 = -pkin(6) + t119;
t173 = t112 - t119;
t214 = t173 * t121;
t118 = cos(pkin(7));
t101 = t118 * pkin(2) + pkin(1);
t203 = rSges(4,2) * t104;
t205 = rSges(4,1) * t105;
t150 = -t203 + t205;
t133 = -t101 - t150;
t49 = (rSges(4,3) - t119) * t121 + t133 * t120;
t111 = t120 * rSges(4,3);
t188 = t121 * t205 + t111;
t50 = -t120 * t119 + (t101 - t203) * t121 + t188;
t213 = t120 * t50 + t121 * t49;
t85 = Icges(5,2) * t100 + t201;
t86 = Icges(5,1) * t99 + t185;
t151 = t100 * t86 - t85 * t99;
t212 = qJD(1) * t151 - t141 * t114;
t136 = Icges(4,5) * t105 - Icges(4,6) * t104;
t68 = -Icges(4,3) * t121 + t120 * t136;
t135 = rSges(3,1) * t118 - rSges(3,2) * sin(pkin(7)) + pkin(1);
t191 = rSges(3,3) + qJ(2);
t59 = t120 * t191 + t121 * t135;
t211 = 2 * m(4);
t210 = 2 * m(5);
t115 = t120 ^ 2;
t116 = t121 ^ 2;
t92 = rSges(4,1) * t104 + rSges(4,2) * t105;
t209 = m(4) * t92;
t208 = t120 / 0.2e1;
t207 = -t121 / 0.2e1;
t66 = -rSges(5,3) * t121 + t120 * t154;
t32 = t120 * t66 + t121 * t67;
t200 = t104 * t72;
t199 = t104 * t73;
t198 = t105 * t70;
t197 = t105 * t71;
t196 = t114 * t85;
t195 = t114 * t86;
t194 = t114 * t99;
t190 = rSges(5,3) - t112;
t170 = qJD(1) * t121;
t171 = qJD(1) * t120;
t189 = rSges(5,3) * t170 + t171 * t206;
t61 = Icges(5,3) * t120 + t121 * t141;
t178 = qJD(1) * t61;
t69 = Icges(4,3) * t120 + t121 * t136;
t177 = qJD(1) * t69;
t176 = t100 * t114;
t175 = t114 * t121;
t174 = t120 * t112;
t172 = t115 + t116;
t168 = qJD(3) * t105;
t167 = qJD(3) * t120;
t166 = t120 * (t67 * qJD(1) - t120 * t132) + t121 * (-t175 * t202 + (-t100 * t171 - t175 * t99) * rSges(5,1) + t189) + t66 * t170;
t165 = t121 * t203;
t163 = t104 * t171;
t162 = -pkin(3) * t104 - t87;
t35 = -qJD(1) * t62 - t121 * t196;
t159 = t114 * t65 + t35;
t36 = qJD(1) * t63 - t120 * t196;
t158 = t114 * t64 + t36;
t37 = -qJD(1) * t64 - t121 * t195;
t157 = -t114 * t63 + t37;
t38 = qJD(1) * t65 - t120 * t195;
t156 = t114 * t62 - t38;
t12 = t120 * t153 - t121 * t60;
t13 = -t121 * t61 + t217;
t14 = t120 * t60 + t215;
t15 = t120 * t61 + t121 * t152;
t84 = Icges(5,5) * t99 + Icges(5,6) * t100;
t127 = t114 * t84;
t33 = -t121 * t127 - t219;
t34 = -t120 * t127 + t178;
t155 = -t121 * ((t121 * t34 + (t13 - t215) * qJD(1)) * t121 + (t12 * qJD(1) + (t100 * t37 - t176 * t63 - t194 * t65 - t35 * t99 + t178) * t120 + (-t33 + t158 * t99 + t156 * t100 + (t152 - t60) * qJD(1)) * t121) * t120) + t120 * ((t120 * t33 + (t14 - t217) * qJD(1)) * t120 + (t15 * qJD(1) + (-t100 * t38 + t176 * t62 + t194 * t64 + t36 * t99 - t219) * t121 + (-t34 - t159 * t99 + t157 * t100 + (t153 + t61) * qJD(1)) * t120) * t121) + (-t12 * t121 + t13 * t120) * t171 + (t15 * t120 - t121 * t14) * t170;
t93 = pkin(3) * t105 + t101;
t134 = -t154 - t93;
t45 = t120 * t134 + t121 * t190;
t88 = t121 * t93;
t46 = t67 + t88 - t174;
t145 = t120 * t46 + t121 * t45;
t54 = t162 * t120;
t55 = t162 * t121;
t144 = t120 * t54 + t121 * t55;
t139 = Icges(4,1) * t104 + t186;
t137 = Icges(4,2) * t105 + t187;
t75 = t142 * t114;
t76 = t143 * t114;
t122 = qJD(1) * t84 + (-t75 - t195) * t99 + (t76 - t196) * t100;
t131 = (t100 * t159 - t120 * t212 + t122 * t121 + t157 * t99) * t208 + (t100 * t158 + t122 * t120 + t212 * t121 - t156 * t99) * t207 + (t100 * t62 + t120 * t151 - t121 * t84 + t64 * t99) * t171 / 0.2e1 + (t100 * t63 + t120 * t84 + t121 * t151 + t65 * t99) * t170 / 0.2e1;
t130 = qJD(3) * t92;
t126 = qJD(3) * t139;
t125 = qJD(3) * t137;
t124 = qJD(3) * (-Icges(4,5) * t104 - Icges(4,6) * t105);
t123 = rSges(4,2) * t163 + rSges(4,3) * t170 - t121 * t130;
t58 = -t120 * t135 + t121 * t191;
t108 = qJD(2) * t121;
t107 = qJD(2) * t120;
t98 = t119 * t171;
t83 = t150 * qJD(3);
t79 = t154 * t114;
t78 = -t165 + t188;
t77 = -rSges(4,3) * t121 + t120 * t150;
t52 = -t121 * t101 - t120 * t173 + t88;
t51 = t214 + (-t101 + t93) * t120;
t48 = -qJD(1) * t59 + t108;
t47 = qJD(1) * t58 + t107;
t40 = t120 * t124 + t177;
t39 = -qJD(1) * t68 + t121 * t124;
t29 = -t87 * t170 - t120 * t79 + (-t104 * t170 - t105 * t167) * pkin(3);
t28 = t87 * t171 - t121 * t79 + (-t121 * t168 + t163) * pkin(3);
t27 = t108 + t98 + t92 * t167 + (t121 * t133 - t111) * qJD(1);
t26 = t107 + (-t119 * t121 + (-t101 - t205) * t120) * qJD(1) + t123;
t21 = t120 * t69 - t146 * t121;
t20 = t120 * t68 - t216;
t19 = -t121 * t69 - t218;
t18 = -t120 * t147 - t121 * t68;
t17 = t108 + t220 * t120 + (-t120 * t190 + t121 * t134) * qJD(1);
t16 = t107 + (-t93 - t204) * t171 + (-qJD(1) * t112 - t220) * t121 + t189;
t11 = t120 * t51 + t121 * t52 + t32;
t10 = -t171 * t67 + t166;
t3 = t120 * (-t120 * t164 + t98) - t116 * t164 + ((t51 - t214) * t121 + (-t52 - t67 - t174) * t120) * qJD(1) + t166;
t1 = [t86 * t176 + t99 * t76 - t85 * t194 + t100 * t75 + (t16 * t46 + t17 * t45) * t210 + (t26 * t50 + t27 * t49) * t211 + 0.2e1 * m(3) * (t47 * t59 + t48 * t58) + (t140 - t137) * t169 + (t139 + t138) * t168; m(5) * (qJD(1) * t145 + t120 * t17 - t121 * t16) + m(4) * (qJD(1) * t213 + t120 * t27 - t121 * t26) + m(3) * (t120 * t48 - t121 * t47 + (t120 * t59 + t121 * t58) * qJD(1)); 0; (-qJD(3) * t147 + t104 * (qJD(1) * t73 - t120 * t126) + t105 * (qJD(1) * t71 - t120 * t125)) * t207 + (-qJD(3) * t146 + t104 * (-qJD(1) * t72 - t121 * t126) + t105 * (-qJD(1) * t70 - t121 * t125)) * t208 + m(5) * (t16 * t54 + t17 * t55 + t28 * t45 + t29 * t46) + m(4) * ((-t120 * t26 - t121 * t27) * t92 - t213 * t83) + (t115 / 0.2e1 + t116 / 0.2e1) * t136 * qJD(3) + ((-t50 * t209 + t199 / 0.2e1 + t197 / 0.2e1) * t121 + (t200 / 0.2e1 + t198 / 0.2e1 + t49 * t209) * t120) * qJD(1) + t131; m(5) * (qJD(1) * t144 + t28 * t120 - t121 * t29); ((t120 * t77 + t121 * t78) * ((qJD(1) * t77 + t123) * t121 + (-t120 * t130 + (-t165 - t78 + t111) * qJD(1)) * t120) + t172 * t92 * t83) * t211 + (t21 * t120 - t121 * t20) * t170 + t120 * ((t120 * t39 + (t20 + t218) * qJD(1)) * t120 + (t21 * qJD(1) + (t168 * t70 + t169 * t72) * t121 + (-t40 + (-t197 - t199) * qJD(3) + (-t147 + t69) * qJD(1)) * t120) * t121) + (t19 * t120 - t121 * t18) * t171 - t121 * ((t121 * t40 + (t19 + t216) * qJD(1)) * t121 + (t18 * qJD(1) + (-t168 * t71 - t169 * t73 + t177) * t120 + (-t39 + (t198 + t200) * qJD(3) - t146 * qJD(1)) * t121) * t120) + (t11 * t3 + t28 * t55 + t29 * t54) * t210 + t155; m(5) * (-t145 * t79 + (-t120 * t16 - t121 * t17 + (t120 * t45 - t121 * t46) * qJD(1)) * t87) + t131; 0; m(5) * (t10 * t11 + t32 * t3 - t144 * t79 + (-t120 * t29 - t121 * t28 + (t120 * t55 - t121 * t54) * qJD(1)) * t87) + t155; (t172 * t79 * t87 + t10 * t32) * t210 + t155;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
