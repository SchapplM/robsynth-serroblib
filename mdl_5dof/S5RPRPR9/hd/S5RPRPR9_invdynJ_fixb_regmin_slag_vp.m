% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:39
% EndTime: 2019-12-31 18:24:43
% DurationCPUTime: 1.28s
% Computational Cost: add. (1007->232), mult. (1985->303), div. (0->0), fcn. (1189->10), ass. (0->143)
t81 = sin(pkin(8));
t64 = t81 * pkin(1) + pkin(6);
t54 = t64 * qJDD(1);
t182 = qJD(2) * qJD(3) + t54;
t135 = qJD(1) * qJD(3);
t84 = sin(qJ(3));
t126 = t84 * t135;
t87 = cos(qJ(3));
t137 = t87 * qJDD(1);
t181 = t126 - t137;
t56 = t64 * qJD(1);
t158 = t87 * qJD(2) - t84 * t56;
t180 = -qJD(4) + t158;
t19 = -qJD(3) * pkin(3) - t180;
t141 = qJD(3) * qJ(4);
t29 = t84 * qJD(2) + t87 * t56;
t21 = -t141 - t29;
t83 = sin(qJ(5));
t145 = t83 * qJD(3);
t153 = qJD(1) * t87;
t86 = cos(qJ(5));
t48 = t86 * t153 + t145;
t152 = qJD(3) * t48;
t125 = t87 * t135;
t138 = t84 * qJDD(1);
t106 = t125 + t138;
t43 = qJDD(5) + t106;
t35 = t86 * t43;
t179 = t152 - t35;
t78 = qJ(1) + pkin(8);
t71 = sin(t78);
t72 = cos(t78);
t178 = g(1) * t72 + g(2) * t71;
t144 = t84 * qJD(1);
t142 = pkin(4) * t144 - t180;
t150 = qJD(3) * t87;
t177 = -t56 * t150 - t182 * t84;
t68 = pkin(4) * t153;
t14 = -t21 + t68;
t63 = qJD(5) + t144;
t89 = -pkin(3) - pkin(7);
t176 = t89 * t43 + (t14 - t68 - t29) * t63;
t148 = qJD(5) * t87;
t103 = t84 * t145 - t86 * t148;
t164 = t83 * t87;
t175 = -t103 * t63 + t43 * t164;
t128 = t83 * t153;
t12 = -qJD(5) * t128 + t83 * qJDD(3) + (qJD(3) * qJD(5) - t181) * t86;
t172 = g(3) * t84;
t77 = g(3) * t87;
t82 = cos(pkin(8));
t171 = t82 * pkin(1);
t170 = pkin(4) + t64;
t11 = -t48 * qJD(5) + t86 * qJDD(3) + t181 * t83;
t169 = t11 * t86;
t168 = t48 * t63;
t50 = t86 * qJD(3) - t128;
t167 = t50 * t63;
t166 = t83 * t43;
t165 = t83 * t84;
t163 = t84 * t12;
t162 = t84 * t86;
t161 = t86 * t87;
t160 = t11 * t84 + t50 * t150;
t149 = qJD(5) * t63;
t129 = t83 * t149;
t132 = t63 * t162;
t159 = qJD(3) * t132 + t87 * t129;
t79 = t84 ^ 2;
t80 = t87 ^ 2;
t157 = t79 - t80;
t156 = qJ(4) * t87;
t155 = t84 * qJ(4);
t117 = t87 * pkin(3) + t155;
t111 = pkin(2) + t117;
t36 = -t111 - t171;
t22 = qJD(1) * t36;
t70 = pkin(3) * t144;
t51 = -qJ(4) * t153 + t70;
t154 = qJD(1) * t51;
t65 = -pkin(2) - t171;
t57 = qJD(1) * t65;
t151 = qJD(3) * t84;
t147 = qJDD(3) * pkin(3);
t146 = t29 * qJD(3);
t143 = t84 * qJD(4);
t140 = qJDD(1) * t36;
t139 = qJDD(3) * t64;
t136 = t87 * qJDD(2);
t133 = qJDD(3) * qJ(4);
t91 = qJD(1) ^ 2;
t130 = t84 * t91 * t87;
t110 = t136 + t177;
t104 = qJDD(4) - t110;
t4 = t106 * pkin(4) + t89 * qJDD(3) + t104;
t27 = t89 * t87 - t155 + t65;
t62 = pkin(3) * t126;
t115 = pkin(7) * t84 - t156;
t98 = t115 * qJD(3) - t143;
t6 = t98 * qJD(1) + t27 * qJDD(1) + t62;
t127 = t86 * t4 - t83 * t6;
t38 = t170 * t87;
t13 = t89 * qJD(3) + t142;
t124 = qJD(5) * t13 + t6;
t15 = t27 * qJD(1);
t123 = qJD(5) * t15 - t4;
t122 = t84 * qJDD(2) - t56 * t151 + t182 * t87;
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t119 = g(1) * t85 - g(2) * t88;
t116 = pkin(3) * t84 - t156;
t2 = t83 * t13 + t86 * t15;
t112 = t63 * t83;
t109 = t119 * pkin(1);
t107 = -t87 * t141 - t143;
t90 = qJD(3) ^ 2;
t105 = g(1) * t71 - g(2) * t72 - t64 * t90;
t102 = -qJD(1) * t57 + t178;
t101 = 0.2e1 * t57 * qJD(3) - t139;
t100 = -0.2e1 * t22 * qJD(3) + t139;
t99 = qJD(3) * qJD(4) + t122 + t133;
t97 = -0.2e1 * qJDD(1) * t65 + t105;
t8 = t104 - t147;
t96 = t99 * t87 + t8 * t84 + (t19 * t87 + t21 * t84) * qJD(3);
t95 = t22 * t144 - t178 * t84 + qJDD(4) - t177 + t77;
t69 = pkin(3) * t151;
t34 = t107 + t69;
t9 = t107 * qJD(1) + t140 + t62;
t94 = -qJD(1) * t34 + t105 - t140 - t9;
t5 = -t181 * pkin(4) + t99;
t93 = -t172 + t5 - t178 * t87 + (t115 * qJD(1) - qJD(5) * t89 + t70) * t63;
t53 = qJDD(3) * t87 - t90 * t84;
t52 = qJDD(3) * t84 + t90 * t87;
t37 = t170 * t84;
t33 = qJD(3) * t38;
t32 = t170 * t151;
t26 = -t71 * t165 + t72 * t86;
t25 = t71 * t162 + t72 * t83;
t24 = t72 * t165 + t71 * t86;
t23 = t72 * t162 - t71 * t83;
t20 = t69 + t98;
t1 = t86 * t13 - t83 * t15;
t3 = [qJDD(1), t119, g(1) * t88 + g(2) * t85, (t81 ^ 2 + t82 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t109, t79 * qJDD(1) + 0.2e1 * t84 * t125, -0.2e1 * t157 * t135 + 0.2e1 * t84 * t137, t52, t53, 0, t101 * t84 + t87 * t97, t101 * t87 - t84 * t97, (t79 + t80) * t54 + t96 - t178, t100 * t84 - t87 * t94, t100 * t87 + t84 * t94, t22 * t34 + t9 * t36 + t109 + (-g(1) * pkin(6) - g(2) * t111) * t72 + (-g(2) * pkin(6) + g(1) * t111) * t71 + t96 * t64, t103 * t50 - t11 * t164, (-t48 * t83 + t50 * t86) * t151 + (-t169 + t12 * t83 + (t48 * t86 + t50 * t83) * qJD(5)) * t87, t160 - t175, -t163 + (-t152 - t35) * t87 + t159, t63 * t150 + t43 * t84, (-t83 * t20 + t86 * t33) * t63 + (-t83 * t27 + t86 * t37) * t43 + t127 * t84 - t32 * t48 + t38 * t12 + t5 * t161 - g(1) * t26 - g(2) * t24 + (t1 * t87 - t14 * t162) * qJD(3) + ((-t86 * t27 - t83 * t37) * t63 - t2 * t84 - t14 * t164) * qJD(5), -t2 * t150 + g(1) * t25 - g(2) * t23 + t38 * t11 - t32 * t50 + (-(qJD(5) * t37 + t20) * t63 - t27 * t43 - t124 * t84 - t14 * t148) * t86 + (-(-qJD(5) * t27 + t33) * t63 - t37 * t43 - t5 * t87 + (t14 * qJD(3) + t123) * t84) * t83; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t53, -t52, 0, -t53, t52, t99 * t84 - t8 * t87 - g(3) + (t19 * t84 - t21 * t87) * qJD(3), 0, 0, 0, 0, 0, t179 * t87 + t159 + t163, t160 + t175; 0, 0, 0, 0, -t130, t157 * t91, t138, t137, qJDD(3), t102 * t84 + t110 + t146 - t77, qJD(3) * t158 + t102 * t87 - t122 + t172, -t116 * qJDD(1), -0.2e1 * t147 - t146 + (-qJDD(2) - t154) * t87 + t95, 0.2e1 * t133 + (-g(3) + t154) * t84 + (0.2e1 * qJD(4) - t158) * qJD(3) + (qJD(1) * t22 - t178) * t87 + t122, -t8 * pkin(3) - g(3) * t117 + t99 * qJ(4) + t178 * t116 + t180 * t21 - t19 * t29 - t22 * t51, -t112 * t50 + t169, (-t12 - t167) * t86 + (-t11 + t168) * t83, -t129 + t35 + (-t63 * t165 - t50 * t87) * qJD(1), -t86 * t149 - t166 + (t48 * t87 - t132) * qJD(1), -t63 * t153, qJ(4) * t12 - t1 * t153 + t142 * t48 + t176 * t86 + t93 * t83, qJ(4) * t11 + t142 * t50 + t2 * t153 - t176 * t83 + t93 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, qJDD(3) + t130, -t79 * t91 - t90, t21 * qJD(3) - t136 - t147 + t95, 0, 0, 0, 0, 0, -t63 * t112 - t179, -t63 ^ 2 * t86 - qJD(3) * t50 - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t48, -t48 ^ 2 + t50 ^ 2, t11 + t168, -t12 + t167, t43, -g(1) * t23 - g(2) * t25 + g(3) * t161 - t14 * t50 + t127 + (-qJD(5) + t63) * t2, g(1) * t24 - g(2) * t26 + t1 * t63 + t14 * t48 - t124 * t86 + (t123 - t77) * t83;];
tau_reg = t3;
