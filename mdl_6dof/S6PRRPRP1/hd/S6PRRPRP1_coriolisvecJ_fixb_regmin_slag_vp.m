% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [6x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:34
% EndTime: 2019-03-08 21:26:41
% DurationCPUTime: 2.28s
% Computational Cost: add. (3204->289), mult. (8274->420), div. (0->0), fcn. (6283->10), ass. (0->165)
t127 = sin(qJ(2));
t122 = sin(pkin(6));
t179 = qJD(1) * t122;
t165 = t127 * t179;
t126 = sin(qJ(3));
t173 = t126 * qJD(3);
t219 = pkin(3) * t173 - t165;
t121 = sin(pkin(11));
t123 = cos(pkin(11));
t129 = cos(qJ(3));
t203 = -qJ(4) - pkin(8);
t157 = qJD(3) * t203;
t91 = qJD(4) * t129 + t126 * t157;
t92 = -qJD(4) * t126 + t129 * t157;
t49 = t121 * t92 + t123 * t91;
t184 = t123 * t129;
t102 = t121 * t126 - t184;
t130 = cos(qJ(2));
t164 = t130 * t179;
t74 = t102 * t164;
t199 = t49 + t74;
t103 = t121 * t129 + t123 * t126;
t95 = t103 * qJD(3);
t98 = t102 * qJD(3);
t218 = pkin(4) * t95 + pkin(9) * t98 + t219;
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t217 = -t125 * t74 + t128 * t218;
t177 = qJD(2) * t126;
t94 = qJD(2) * t184 - t121 * t177;
t90 = qJD(5) - t94;
t155 = t125 * t90;
t96 = t103 * qJD(2);
t80 = qJD(3) * t125 + t128 * t96;
t216 = t80 * t155;
t175 = qJD(5) * t125;
t88 = qJD(2) * t95;
t215 = t128 * t88 - t90 * t175;
t174 = qJD(5) * t128;
t166 = -pkin(3) * t129 - pkin(2);
t61 = pkin(4) * t102 - pkin(9) * t103 + t166;
t214 = t218 * t125 + t128 * t199 + t61 * t174;
t105 = qJD(2) * pkin(8) + t165;
t146 = qJD(4) + t164;
t124 = cos(pkin(6));
t178 = qJD(1) * t124;
t163 = t126 * t178;
t213 = (-t105 * t129 - t163) * qJD(3) + (-qJ(4) * qJD(3) * t129 - t126 * t146) * qJD(2);
t153 = qJ(4) * qJD(2) + t105;
t72 = t129 * t153 + t163;
t195 = t123 * t72;
t112 = t129 * t178;
t71 = -t126 * t153 + t112;
t66 = qJD(3) * pkin(3) + t71;
t25 = t121 * t66 + t195;
t22 = qJD(3) * pkin(9) + t25;
t89 = t166 * qJD(2) + qJD(4) - t164;
t36 = -pkin(4) * t94 - pkin(9) * t96 + t89;
t14 = t125 * t36 + t128 * t22;
t40 = (-t105 * t126 + t112) * qJD(3) + (-qJ(4) * t173 + t129 * t146) * qJD(2);
t17 = t213 * t121 + t123 * t40;
t171 = qJD(2) * qJD(3);
t158 = t129 * t171;
t159 = t126 * t171;
t143 = -t121 * t159 + t123 * t158;
t176 = qJD(2) * t127;
t162 = t122 * t176;
t93 = pkin(3) * t159 + qJD(1) * t162;
t33 = t88 * pkin(4) - pkin(9) * t143 + t93;
t31 = t128 * t33;
t134 = -qJD(5) * t14 - t125 * t17 + t31;
t172 = t128 * qJD(3);
t41 = -qJD(5) * t172 - t128 * t143 + t96 * t175;
t1 = pkin(5) * t88 + t41 * qJ(6) - t80 * qJD(6) + t134;
t78 = t125 * t96 - t172;
t8 = -t78 * qJ(6) + t14;
t212 = t8 * t90 + t1;
t211 = t80 ^ 2;
t13 = -t125 * t22 + t128 * t36;
t7 = -t80 * qJ(6) + t13;
t5 = pkin(5) * t90 + t7;
t210 = t5 - t7;
t145 = qJ(6) * t98 - qJD(6) * t103;
t107 = t203 * t126;
t108 = t203 * t129;
t76 = t107 * t121 - t108 * t123;
t67 = t128 * t76;
t209 = pkin(5) * t95 - t125 * t49 + t145 * t128 + (-t67 + (qJ(6) * t103 - t61) * t125) * qJD(5) + t217;
t160 = t103 * t174;
t208 = -qJ(6) * t160 + (-qJD(5) * t76 + t145) * t125 + t214;
t115 = pkin(3) * t121 + pkin(9);
t181 = qJ(6) + t115;
t154 = qJD(5) * t181;
t189 = qJ(6) * t128;
t62 = t121 * t72;
t29 = t123 * t71 - t62;
t50 = pkin(3) * t177 + pkin(4) * t96 - pkin(9) * t94;
t46 = t128 * t50;
t207 = -pkin(5) * t96 - t128 * t154 + t94 * t189 - t46 + (-qJD(6) + t29) * t125;
t206 = t78 * t94;
t205 = t78 * t96;
t204 = t80 * t96;
t202 = t125 * t50 + t128 * t29;
t137 = t125 * t143;
t42 = qJD(5) * t80 + t137;
t201 = -t125 * t42 - t78 * t174;
t200 = -t103 * t164 + t121 * t91 - t123 * t92;
t198 = t125 * t61 + t67;
t192 = t125 * t94;
t197 = qJ(6) * t192 + qJD(6) * t128 - t125 * t154 - t202;
t196 = qJD(2) * pkin(2);
t194 = t125 * t41;
t193 = t125 * t88;
t191 = t125 * t98;
t190 = t128 * t98;
t188 = t103 * t125;
t187 = t122 * t127;
t186 = t122 * t130;
t132 = qJD(2) ^ 2;
t185 = t122 * t132;
t131 = qJD(3) ^ 2;
t183 = t131 * t126;
t182 = t131 * t129;
t180 = t126 ^ 2 - t129 ^ 2;
t167 = t127 * t185;
t116 = -pkin(3) * t123 - pkin(4);
t161 = qJD(2) * t186;
t16 = t121 * t40 - t123 * t213;
t27 = t121 * t71 + t195;
t24 = t123 * t66 - t62;
t156 = t128 * t90;
t75 = -t123 * t107 - t108 * t121;
t152 = t126 * t161;
t151 = t129 * t161;
t150 = qJD(5) * t115 * t90 + t16;
t140 = t125 * t33 + t128 * t17 + t36 * t174 - t22 * t175;
t2 = -t42 * qJ(6) - t78 * qJD(6) + t140;
t149 = -t90 * t5 + t2;
t147 = t16 * t103 - t76 * t88;
t144 = t90 * t192 + t215;
t6 = t42 * pkin(5) + t16;
t21 = -qJD(3) * pkin(4) - t24;
t141 = t124 * t129 - t126 * t187;
t99 = t124 * t126 + t129 * t187;
t56 = t121 * t141 + t123 * t99;
t38 = -t125 * t56 - t128 * t186;
t142 = t125 * t186 - t128 * t56;
t139 = t196 * qJD(2);
t138 = -t115 * t88 + t90 * t21;
t135 = -0.2e1 * qJD(3) * t196;
t101 = t181 * t128;
t100 = t181 * t125;
t77 = t78 ^ 2;
t70 = -qJD(3) * t99 - t152;
t69 = qJD(3) * t141 + t151;
t58 = t128 * t61;
t55 = t121 * t99 - t123 * t141;
t28 = t121 * t70 + t123 * t69;
t26 = t121 * t69 - t123 * t70;
t20 = -qJ(6) * t188 + t198;
t19 = t78 * pkin(5) + qJD(6) + t21;
t18 = pkin(5) * t102 - t103 * t189 - t125 * t76 + t58;
t12 = qJD(5) * t142 - t125 * t28 + t128 * t162;
t11 = qJD(5) * t38 + t125 * t162 + t128 * t28;
t3 = [0, 0, -t167, -t130 * t185, 0, 0, 0, 0, 0, -t129 * t167 + (t70 - t152) * qJD(3), t126 * t167 + (-t69 - t151) * qJD(3), t143 * t55 + t26 * t96 + t28 * t94 - t56 * t88, t16 * t55 + t17 * t56 - t24 * t26 + t25 * t28 + (-t130 * t93 + t176 * t89) * t122, 0, 0, 0, 0, 0, t12 * t90 + t26 * t78 + t38 * t88 + t55 * t42, -t11 * t90 + t142 * t88 + t26 * t80 - t55 * t41, -t11 * t78 - t12 * t80 + t142 * t42 + t38 * t41, t1 * t38 + t11 * t8 + t12 * t5 - t142 * t2 + t19 * t26 + t55 * t6; 0, 0, 0, 0, 0.2e1 * t126 * t158, -0.2e1 * t180 * t171, t182, -t183, 0, -pkin(8) * t182 + t126 * t135, pkin(8) * t183 + t129 * t135, -t17 * t102 + t75 * t143 + t199 * t94 + t200 * t96 + t24 * t98 - t25 * t95 + t147, t16 * t75 + t93 * t166 + t17 * t76 + t199 * t25 - t200 * t24 + t219 * t89, -t80 * t190 + (-t128 * t41 - t175 * t80) * t103 -(-t125 * t80 - t128 * t78) * t98 + (t194 - t128 * t42 + (t125 * t78 - t128 * t80) * qJD(5)) * t103, -t102 * t41 + t215 * t103 - t90 * t190 + t80 * t95, t90 * t191 - t102 * t42 - t78 * t95 + (-t174 * t90 - t193) * t103, t102 * t88 + t90 * t95, t58 * t88 + (-t174 * t22 + t31) * t102 + t13 * t95 + t75 * t42 + t21 * t160 + (-t174 * t76 + t217) * t90 + t200 * t78 + ((-qJD(5) * t61 - t49) * t90 + (-qJD(5) * t36 - t17) * t102 - t21 * t98 + t147) * t125, -t198 * t88 - t140 * t102 - t14 * t95 - t75 * t41 - t21 * t190 + (t175 * t76 - t214) * t90 + t200 * t80 + (t16 * t128 - t175 * t21) * t103, t18 * t41 - t20 * t42 - (-t125 * t8 - t128 * t5) * t98 - t209 * t80 - t208 * t78 + (-t1 * t128 - t125 * t2 + (t125 * t5 - t128 * t8) * qJD(5)) * t103, t2 * t20 + t1 * t18 + t6 * (pkin(5) * t188 + t75) + t208 * t8 + t209 * t5 + ((t160 - t191) * pkin(5) + t200) * t19; 0, 0, 0, 0, -t126 * t132 * t129, t180 * t132, 0, 0, 0, t126 * t139, t129 * t139 (t25 - t27) * t96 + (-t29 + t24) * t94 + (-t121 * t88 - t123 * t143) * pkin(3), t24 * t27 - t25 * t29 + (t121 * t17 - t123 * t16 - t177 * t89) * pkin(3), t156 * t80 - t194 (-t41 + t206) * t128 - t216 + t201, t156 * t90 + t193 - t204, t144 + t205, -t90 * t96, t116 * t42 - t13 * t96 - t27 * t78 - t46 * t90 - t150 * t128 + (t29 * t90 + t138) * t125, -t116 * t41 + t150 * t125 + t138 * t128 + t14 * t96 + t202 * t90 - t27 * t80, -t100 * t41 - t101 * t42 - t212 * t125 + t149 * t128 - t197 * t78 - t207 * t80, t2 * t101 - t1 * t100 + t6 * (-pkin(5) * t128 + t116) + t197 * t8 + t207 * t5 + (pkin(5) * t155 - t27) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 ^ 2 - t96 ^ 2, t24 * t96 - t25 * t94 + t93, 0, 0, 0, 0, 0, t144 - t205, -t128 * t90 ^ 2 - t193 - t204 (t41 + t206) * t128 + t216 + t201, t149 * t125 + t212 * t128 - t19 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t78, -t77 + t211, t78 * t90 - t41, -t137 + (-qJD(5) + t90) * t80, t88, t14 * t90 - t21 * t80 + t134, t13 * t90 + t21 * t78 - t140, pkin(5) * t41 - t210 * t78, t210 * t8 + (-t19 * t80 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 - t211, t5 * t80 + t8 * t78 + t6;];
tauc_reg  = t3;
