% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:15
% EndTime: 2019-12-31 19:51:19
% DurationCPUTime: 1.26s
% Computational Cost: add. (2242->245), mult. (3267->280), div. (0->0), fcn. (2162->12), ass. (0->153)
t121 = sin(pkin(8));
t122 = cos(pkin(8));
t184 = t121 ^ 2 + t122 ^ 2;
t115 = qJDD(1) + qJDD(2);
t119 = qJD(1) + qJD(2);
t125 = sin(qJ(2));
t179 = qJDD(1) * t125;
t127 = cos(qJ(2));
t181 = qJD(2) * t127;
t54 = t115 * qJ(3) + t119 * qJD(3) + (qJD(1) * t181 + t179) * pkin(1);
t219 = t184 * t54;
t120 = qJ(1) + qJ(2);
t111 = cos(t120);
t103 = g(2) * t111;
t205 = t125 * pkin(1);
t173 = qJD(1) * t205;
t203 = t127 * pkin(1);
t195 = qJD(2) * t173 - qJDD(1) * t203;
t166 = qJDD(3) + t195;
t206 = t115 * pkin(2);
t62 = t166 - t206;
t218 = t62 + t103;
t207 = cos(qJ(4));
t169 = t207 * t122;
t124 = sin(qJ(4));
t189 = t124 * t121;
t140 = t169 - t189;
t72 = t207 * t121 + t124 * t122;
t183 = qJD(1) * t127;
t149 = -pkin(1) * t183 + qJD(3);
t100 = t122 * pkin(3) + pkin(2);
t61 = -t100 * t119 + t149;
t171 = t119 * t189;
t63 = -t119 * t169 + t171;
t65 = t72 * t119;
t18 = t63 * pkin(4) - t65 * qJ(5) + t61;
t161 = qJD(4) * t207;
t153 = t122 * t161;
t180 = qJD(4) * t124;
t167 = t121 * t180;
t67 = -t153 + t167;
t175 = t72 * t115 + t119 * t153;
t28 = t119 * t167 - t175;
t151 = t140 * t115;
t68 = t72 * qJD(4);
t29 = t119 * t68 - t151;
t46 = -t100 * t115 + t166;
t130 = t29 * pkin(4) + t28 * qJ(5) + t46;
t8 = -t65 * qJD(5) + t130;
t217 = t18 * t67 - t8 * t72;
t216 = -t140 * t8 + t18 * t68;
t215 = -t140 * t46 + t61 * t68;
t214 = t46 * t72 - t61 * t67;
t213 = t119 * t184;
t110 = sin(t120);
t212 = g(1) * t111 + g(2) * t110;
t104 = g(1) * t110;
t211 = t103 - t104;
t118 = pkin(8) + qJ(4);
t108 = sin(t118);
t109 = cos(t118);
t164 = pkin(7) * t115 + t54;
t37 = t164 * t121;
t38 = t164 * t122;
t74 = t119 * qJ(3) + t173;
t163 = pkin(7) * t119 + t74;
t55 = t163 * t121;
t177 = -t124 * t37 - t55 * t161 + t207 * t38;
t210 = -g(3) * t108 - t212 * t109 + t177;
t192 = t108 * t111;
t193 = t108 * t110;
t198 = g(1) * t193 - g(2) * t192;
t123 = -pkin(7) - qJ(3);
t85 = t123 * t121;
t112 = t122 * pkin(7);
t86 = t122 * qJ(3) + t112;
t141 = -t124 * t86 + t207 * t85;
t201 = t141 * qJD(4) + t149 * t140;
t53 = t124 * t85 + t207 * t86;
t209 = t201 * qJD(4) + t53 * qJDD(4) + t198;
t208 = t65 ^ 2;
t126 = sin(qJ(1));
t204 = t126 * pkin(1);
t202 = t65 * t63;
t200 = t53 * qJD(4) + t149 * t72;
t199 = t218 * t121;
t56 = t163 * t122;
t197 = t124 * t56;
t196 = t111 * pkin(2) + t110 * qJ(3);
t194 = qJDD(4) * pkin(4);
t191 = t109 * t111;
t190 = t111 * t123;
t23 = -t124 * t55 + t207 * t56;
t187 = t23 * qJD(4);
t22 = -t207 * t55 - t197;
t186 = qJD(5) - t22;
t182 = qJD(2) * t125;
t178 = qJDD(4) * qJ(5);
t176 = pkin(4) * t191 + qJ(5) * t192 + t111 * t100;
t174 = pkin(1) * t182;
t168 = t119 * t182;
t165 = -t110 * pkin(2) + t111 * qJ(3);
t160 = t22 + t197;
t159 = t124 * t38 + t56 * t161 - t55 * t180 + t207 * t37;
t158 = t184 * t115;
t157 = -t212 + t219;
t156 = t119 * t173;
t155 = t195 + t211;
t154 = -g(2) * t191 + t109 * t104;
t24 = t68 * pkin(4) + t67 * qJ(5) - t72 * qJD(5);
t152 = -t24 + t173;
t146 = t109 * pkin(4) + t108 * qJ(5);
t19 = -qJD(4) * pkin(4) + t186;
t20 = qJD(4) * qJ(5) + t23;
t6 = t178 + (qJD(5) - t197) * qJD(4) + t177;
t7 = qJDD(5) + t159 - t194;
t143 = t140 * t6 - t19 * t67 - t20 * t68 + t7 * t72 - t212;
t99 = qJ(3) + t205;
t69 = (-pkin(7) - t99) * t121;
t70 = t122 * t99 + t112;
t142 = -t124 * t70 + t207 * t69;
t42 = t124 * t69 + t207 * t70;
t139 = -t156 - t206;
t106 = -pkin(2) - t203;
t138 = pkin(1) * t168 + t106 * t115;
t93 = pkin(1) * t181 + qJD(3);
t14 = t142 * qJD(4) + t140 * t93;
t137 = t14 * qJD(4) + t42 * qJDD(4) + t198;
t136 = g(1) * t192 + g(2) * t193 - g(3) * t109 - t159;
t43 = -pkin(4) * t140 - t72 * qJ(5) - t100;
t15 = t42 * qJD(4) + t72 * t93;
t135 = -t15 * qJD(4) + qJDD(4) * t142 + t154;
t134 = t149 * t184;
t133 = -t200 * qJD(4) + qJDD(4) * t141 + t154;
t132 = t18 * t65 + qJDD(5) - t136;
t131 = (-g(1) * (-t100 - t146) + g(2) * t123) * t110;
t129 = 0.2e1 * t65 * qJD(4) - t151;
t128 = cos(qJ(1));
t113 = t128 * pkin(1);
t92 = t122 * t104;
t84 = -t100 - t203;
t73 = -t119 * pkin(2) + t149;
t60 = t63 ^ 2;
t45 = -t68 * qJD(4) + qJDD(4) * t140;
t44 = -t67 * qJD(4) + t72 * qJDD(4);
t36 = t43 - t203;
t27 = t65 * pkin(4) + t63 * qJ(5);
t21 = t24 + t174;
t17 = (t63 - t171) * qJD(4) + t175;
t16 = (t63 + t171) * qJD(4) - t175;
t9 = -t28 * t72 - t65 * t67;
t1 = -t140 * t28 - t72 * t29 + t67 * t63 - t65 * t68;
t2 = [qJDD(1), g(1) * t126 - g(2) * t128, g(1) * t128 + g(2) * t126, t115, (t115 * t127 - t168) * pkin(1) - t155, ((-qJDD(1) - t115) * t125 + (-qJD(1) - t119) * t181) * pkin(1) + t212, t92 + (-t138 - t218) * t122, (t138 - t104) * t121 + t199, t99 * t158 + t93 * t213 + t157, t62 * t106 + t73 * t174 - g(1) * (t165 - t204) - g(2) * (t113 + t196) + t184 * (t54 * t99 + t74 * t93), t9, t1, t44, t45, 0, t63 * t174 + t84 * t29 + t135 + t215, t65 * t174 - t84 * t28 - t137 + t214, t21 * t63 + t36 * t29 + t135 + t216, -t14 * t63 + t142 * t28 + t15 * t65 - t42 * t29 + t143, -t21 * t65 + t36 * t28 + t137 + t217, t6 * t42 + t20 * t14 + t8 * t36 + t18 * t21 - t7 * t142 + t19 * t15 - g(1) * (-t190 - t204) - g(2) * (t113 + t176) + t131; 0, 0, 0, t115, -t155 + t156, (-t179 + (-qJD(2) + t119) * t183) * pkin(1) + t212, t92 + (-t139 - t218) * t122, (t139 - t104) * t121 + t199, qJ(3) * t158 + t134 * t119 + t157, -t62 * pkin(2) - g(1) * t165 - g(2) * t196 + qJ(3) * t219 + t134 * t74 - t73 * t173, t9, t1, t44, t45, 0, -t100 * t29 - t63 * t173 + t133 + t215, t100 * t28 - t65 * t173 - t209 + t214, -t152 * t63 + t43 * t29 + t133 + t216, t141 * t28 + t200 * t65 - t201 * t63 - t53 * t29 + t143, t152 * t65 + t43 * t28 + t209 + t217, g(1) * t190 - g(2) * t176 - t141 * t7 - t152 * t18 + t200 * t19 + t201 * t20 + t8 * t43 + t6 * t53 + t131; 0, 0, 0, 0, 0, 0, -t122 * t115, t121 * t115, -t184 * t119 ^ 2, -t74 * t213 + t211 + t62, 0, 0, 0, 0, 0, t129, -t16, t129, -t60 - t208, t16, t20 * t63 + (-qJD(5) - t19) * t65 + t130 + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t60 + t208, t17, t151, qJDD(4), -t61 * t65 + t136 + t187, t160 * qJD(4) + t61 * t63 - t210, -t27 * t63 - t132 + t187 + 0.2e1 * t194, pkin(4) * t28 - t29 * qJ(5) + (t20 - t23) * t65 + (t19 - t186) * t63, 0.2e1 * t178 - t18 * t63 + t27 * t65 + (0.2e1 * qJD(5) - t160) * qJD(4) + t210, -t7 * pkin(4) - g(3) * t146 + t6 * qJ(5) - t18 * t27 + t186 * t20 - t19 * t23 + t212 * (pkin(4) * t108 - qJ(5) * t109); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t202, t17, -qJD(4) ^ 2 - t208, -t20 * qJD(4) + t132 - t194;];
tau_reg = t2;
