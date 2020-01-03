% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP7
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
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:34
% EndTime: 2019-12-31 20:01:41
% DurationCPUTime: 2.30s
% Computational Cost: add. (3548->367), mult. (8195->466), div. (0->0), fcn. (5791->10), ass. (0->172)
t131 = cos(qJ(2));
t196 = cos(pkin(8));
t168 = t196 * t131;
t107 = qJD(1) * t168;
t125 = sin(pkin(8));
t128 = sin(qJ(2));
t186 = qJD(1) * t128;
t83 = -t125 * t186 + t107;
t75 = qJD(4) - t83;
t122 = qJ(2) + pkin(8);
t116 = sin(t122);
t129 = sin(qJ(1));
t132 = cos(qJ(1));
t162 = g(1) * t132 + g(2) * t129;
t230 = t116 * t162;
t229 = g(1) * t129 - g(2) * t132;
t117 = cos(t122);
t216 = t116 * pkin(7);
t228 = t117 * pkin(3) + t216;
t207 = qJ(3) + pkin(6);
t102 = t207 * t131;
t97 = qJD(1) * t102;
t88 = t125 * t97;
t201 = qJD(2) * pkin(2);
t101 = t207 * t128;
t96 = qJD(1) * t101;
t91 = -t96 + t201;
t58 = t196 * t91 - t88;
t48 = -qJD(2) * pkin(3) - t58;
t127 = sin(qJ(4));
t130 = cos(qJ(4));
t94 = t125 * t131 + t128 * t196;
t85 = t94 * qJD(1);
t68 = -t130 * qJD(2) + t127 * t85;
t70 = t127 * qJD(2) + t130 * t85;
t16 = t68 * pkin(4) - t70 * qJ(5) + t48;
t214 = t125 * pkin(2);
t110 = pkin(7) + t214;
t180 = t128 * qJDD(1);
t156 = qJDD(1) * t168 - t125 * t180;
t84 = t94 * qJD(2);
t55 = qJD(1) * t84 + qJDD(4) - t156;
t200 = t110 * t55;
t227 = t16 * t75 - t200;
t182 = qJD(1) * qJD(2);
t172 = t128 * t182;
t136 = qJD(2) * t107 + qJDD(1) * t94 - t125 * t172;
t226 = t70 ^ 2;
t225 = t75 ^ 2;
t224 = t55 * pkin(4);
t213 = t131 * pkin(2);
t115 = pkin(1) + t213;
t100 = -qJD(1) * t115 + qJD(3);
t34 = -t83 * pkin(3) - t85 * pkin(7) + t100;
t173 = t196 * t97;
t59 = t125 * t91 + t173;
t49 = qJD(2) * pkin(7) + t59;
t15 = t127 * t34 + t130 * t49;
t9 = t75 * qJ(5) + t15;
t223 = t9 * t75;
t219 = g(3) * t116;
t218 = g(3) * t117;
t217 = g(3) * t131;
t212 = t15 * t75;
t211 = t68 * t83;
t210 = t70 * t68;
t169 = t70 * t75;
t209 = t70 * t85;
t208 = t85 * t68;
t183 = qJD(4) * t130;
t135 = -t130 * qJDD(2) + t127 * t136;
t29 = t70 * qJD(4) + t135;
t206 = -t127 * t29 - t68 * t183;
t170 = qJD(2) * t207;
t81 = -t128 * qJD(3) - t131 * t170;
t54 = qJDD(2) * pkin(2) + qJD(1) * t81 - qJDD(1) * t101;
t80 = t131 * qJD(3) - t128 * t170;
t61 = qJD(1) * t80 + qJDD(1) * t102;
t25 = t125 * t54 + t196 * t61;
t41 = pkin(2) * t186 + t85 * pkin(3) - t83 * pkin(7);
t63 = -t196 * t96 - t88;
t205 = t127 * t41 + t130 * t63;
t147 = -t125 * t128 + t168;
t57 = -pkin(3) * t147 - t94 * pkin(7) - t115;
t67 = -t125 * t101 + t102 * t196;
t204 = t127 * t57 + t130 * t67;
t157 = pkin(4) * t127 - qJ(5) * t130;
t62 = -t125 * t96 + t173;
t203 = -t127 * qJD(5) + t75 * t157 - t62;
t189 = t132 * t130;
t192 = t129 * t130;
t202 = (g(1) * t189 + g(2) * t192) * t116;
t45 = t127 * t55;
t199 = t127 * t75;
t46 = t130 * t55;
t181 = qJD(2) * qJD(4);
t184 = qJD(4) * t127;
t28 = -t127 * qJDD(2) + t184 * t85 + (-t136 - t181) * t130;
t198 = t28 * t127;
t197 = t55 * qJ(5);
t195 = qJD(4) * t94;
t194 = qJDD(1) * pkin(1);
t193 = t129 * t127;
t191 = t132 * t207;
t190 = t132 * t127;
t14 = -t127 * t49 + t130 * t34;
t188 = qJD(5) - t14;
t123 = t128 ^ 2;
t187 = -t131 ^ 2 + t123;
t185 = qJD(4) * t110;
t179 = t131 * qJDD(1);
t178 = pkin(2) * t172 + qJDD(3);
t177 = t128 * t201;
t176 = t94 * t184;
t175 = t94 * t183;
t174 = t196 * pkin(2);
t39 = t125 * t80 - t196 * t81;
t60 = -qJD(2) * t85 + t156;
t17 = -pkin(2) * t179 - t60 * pkin(3) - pkin(7) * t136 + t178 - t194;
t23 = qJDD(2) * pkin(7) + t25;
t167 = -t127 * t23 + t130 * t17 - t49 * t183 - t34 * t184;
t66 = t196 * t101 + t125 * t102;
t76 = t117 * t193 + t189;
t78 = t117 * t190 - t192;
t165 = -g(1) * t76 + g(2) * t78;
t77 = t117 * t192 - t190;
t79 = t117 * t189 + t193;
t164 = g(1) * t77 - g(2) * t79;
t24 = -t125 * t61 + t196 * t54;
t22 = -qJDD(2) * pkin(3) - t24;
t3 = t29 * pkin(4) + t28 * qJ(5) - t70 * qJD(5) + t22;
t87 = t147 * qJD(2);
t163 = t16 * t87 + t3 * t94;
t8 = -t75 * pkin(4) + t188;
t160 = -t127 * t9 + t130 * t8;
t159 = t22 * t94 + t48 * t87;
t158 = t55 * t94 + t75 * t87;
t155 = t45 + (-t130 * t83 + t183) * t75;
t154 = -t184 * t75 + t83 * t199 + t46;
t153 = t130 * pkin(4) + t127 * qJ(5) + pkin(3);
t152 = t185 * t75 + t218;
t151 = -0.2e1 * pkin(1) * t182 - pkin(6) * qJDD(2);
t150 = t127 * t17 + t130 * t23 + t34 * t183 - t184 * t49;
t40 = t125 * t81 + t196 * t80;
t42 = t84 * pkin(3) - t87 * pkin(7) + t177;
t149 = t127 * t42 + t130 * t40 + t57 * t183 - t184 * t67;
t148 = t48 * t75 - t200;
t146 = -t152 - t3;
t144 = -qJDD(1) * t115 + t178;
t143 = g(1) * t78 + g(2) * t76 + t127 * t219 + t167;
t133 = qJD(2) ^ 2;
t142 = -pkin(6) * t133 + 0.2e1 * t194 + t229;
t134 = qJD(1) ^ 2;
t141 = pkin(1) * t134 - pkin(6) * qJDD(1) + t162;
t139 = t16 * t70 + qJDD(5) - t143;
t138 = -g(1) * t79 - g(2) * t77 - t130 * t219 + t150;
t111 = -t174 - pkin(3);
t105 = t132 * t115;
t92 = -t174 - t153;
t31 = t70 * pkin(4) + t68 * qJ(5);
t30 = t157 * t94 + t66;
t19 = pkin(4) * t147 + t127 * t67 - t130 * t57;
t18 = -qJ(5) * t147 + t204;
t11 = -t85 * pkin(4) + t127 * t63 - t130 * t41;
t10 = t85 * qJ(5) + t205;
t7 = t68 * t75 - t28;
t6 = (pkin(4) * t87 + qJ(5) * t195) * t127 + (-qJ(5) * t87 + (pkin(4) * qJD(4) - qJD(5)) * t94) * t130 + t39;
t5 = -t84 * pkin(4) + qJD(4) * t204 + t127 * t40 - t130 * t42;
t4 = t84 * qJ(5) - qJD(5) * t147 + t149;
t2 = qJDD(5) - t167 - t224;
t1 = t75 * qJD(5) + t150 + t197;
t12 = [qJDD(1), t229, t162, t123 * qJDD(1) + 0.2e1 * t131 * t172, 0.2e1 * t128 * t179 - 0.2e1 * t182 * t187, qJDD(2) * t128 + t133 * t131, qJDD(2) * t131 - t133 * t128, 0, t128 * t151 + t131 * t142, -t128 * t142 + t131 * t151, t136 * t66 + t147 * t25 - t24 * t94 + t39 * t85 + t40 * t83 - t58 * t87 - t59 * t84 + t67 * t60 - t162, t25 * t67 + t59 * t40 - t24 * t66 - t58 * t39 - t144 * t115 + t100 * t177 - g(1) * (-t129 * t115 + t191) - g(2) * (t129 * t207 + t105), -t70 * t176 + (-t28 * t94 + t70 * t87) * t130, (-t127 * t70 - t130 * t68) * t87 + (t198 - t130 * t29 + (t127 * t68 - t130 * t70) * qJD(4)) * t94, t130 * t158 + t147 * t28 - t176 * t75 + t70 * t84, -t127 * t158 + t147 * t29 - t175 * t75 - t68 * t84, -t147 * t55 + t75 * t84, -t167 * t147 + t14 * t84 + t39 * t68 + t66 * t29 + ((-qJD(4) * t67 + t42) * t75 + t57 * t55 + t48 * t195) * t130 + ((-qJD(4) * t57 - t40) * t75 - t67 * t55 + t159) * t127 + t164, t130 * t159 + t147 * t150 - t149 * t75 - t15 * t84 - t176 * t48 - t204 * t55 - t66 * t28 + t39 * t70 + t165, t127 * t163 + t147 * t2 + t16 * t175 - t19 * t55 + t30 * t29 - t5 * t75 + t6 * t68 - t8 * t84 + t164, -t18 * t29 - t19 * t28 - t4 * t68 + t5 * t70 + t160 * t87 + t229 * t116 + (-t1 * t127 + t2 * t130 + (-t127 * t8 - t130 * t9) * qJD(4)) * t94, -t1 * t147 - t130 * t163 + t16 * t176 + t18 * t55 + t30 * t28 + t4 * t75 - t6 * t70 + t9 * t84 - t165, t1 * t18 + t9 * t4 + t3 * t30 + t16 * t6 + t2 * t19 + t8 * t5 - g(1) * (-t77 * pkin(4) - t76 * qJ(5) + t191) - g(2) * (t79 * pkin(4) + t78 * qJ(5) + t228 * t132 + t105) + (-g(1) * (-t115 - t228) - g(2) * t207) * t129; 0, 0, 0, -t128 * t134 * t131, t187 * t134, t180, t179, qJDD(2), t128 * t141 - t217, g(3) * t128 + t131 * t141, -t136 * t174 + t214 * t60 - (-t59 + t62) * t85 + (-t63 + t58) * t83, t58 * t62 - t59 * t63 + (t196 * t24 - t217 + t125 * t25 + (-qJD(1) * t100 + t162) * t128) * pkin(2), t130 * t169 - t198, (-t28 + t211) * t130 - t70 * t199 + t206, t155 - t209, t154 + t208, -t75 * t85, t111 * t29 - t14 * t85 - t62 * t68 + (-t218 - t22 + (-t41 - t185) * t75) * t130 + (t63 * t75 + t148) * t127 + t202, -t111 * t28 + t205 * t75 + t15 * t85 - t62 * t70 + t148 * t130 + (t152 + t22 - t230) * t127, t11 * t75 + t227 * t127 + t146 * t130 + t203 * t68 + t92 * t29 + t8 * t85 + t202, -t219 + t10 * t68 - t11 * t70 - t162 * t117 + (-t110 * t29 - t8 * t83 + t1 + (t110 * t70 + t8) * qJD(4)) * t130 + (-t110 * t28 + t83 * t9 + t2 + (t110 * t68 - t9) * qJD(4)) * t127, -t10 * t75 + t92 * t28 - t9 * t85 - t203 * t70 - t227 * t130 + (t146 + t230) * t127, t3 * t92 - t9 * t10 - t8 * t11 - g(3) * (t213 + t216) + t203 * t16 - t153 * t218 + (qJD(4) * t160 + t1 * t130 + t2 * t127) * t110 + t162 * (pkin(2) * t128 - pkin(7) * t117 + t116 * t153); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83 ^ 2 - t85 ^ 2, t58 * t85 - t59 * t83 + t144 - t229, 0, 0, 0, 0, 0, t154 - t208, -t130 * t225 - t209 - t45, -t199 * t75 - t208 + t46, (t28 + t211) * t130 + t127 * t169 + t206, t155 + t209, -t16 * t85 + (-t2 + t223) * t130 + (t75 * t8 + t1) * t127 - t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, -t68 ^ 2 + t226, t7, -t127 * t181 - t183 * t85 - t135 + t169, t55, -t48 * t70 + t143 + t212, t14 * t75 + t48 * t68 - t138, -t31 * t68 - t139 + t212 + 0.2e1 * t224, pkin(4) * t28 - t29 * qJ(5) + (-t15 + t9) * t70 + (t8 - t188) * t68, 0.2e1 * t197 - t16 * t68 + t31 * t70 + (0.2e1 * qJD(5) - t14) * t75 + t138, t1 * qJ(5) - t2 * pkin(4) - t16 * t31 - t8 * t15 - g(1) * (-t78 * pkin(4) + t79 * qJ(5)) - g(2) * (-t76 * pkin(4) + t77 * qJ(5)) + t188 * t9 + t157 * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t60 + t210, t7, -t225 - t226, t139 - t223 - t224;];
tau_reg = t12;
