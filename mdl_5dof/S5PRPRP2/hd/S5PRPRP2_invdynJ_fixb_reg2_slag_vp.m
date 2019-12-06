% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:12
% EndTime: 2019-12-05 15:31:15
% DurationCPUTime: 1.37s
% Computational Cost: add. (1201->247), mult. (2530->323), div. (0->0), fcn. (1655->6), ass. (0->154)
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t56 = -t87 * pkin(3) - t86 * pkin(6) - pkin(2);
t90 = cos(qJ(4));
t153 = t87 * t90;
t66 = qJ(3) * t153;
t89 = sin(qJ(4));
t27 = t89 * t56 + t66;
t172 = qJD(4) * t27;
t125 = qJD(2) * qJD(4);
t112 = t90 * t125;
t95 = qJDD(2) * t89 + t112;
t171 = t95 * t86;
t126 = qJD(2) * qJD(3);
t128 = qJ(3) * qJDD(2);
t97 = t126 + t128;
t154 = t87 * t89;
t82 = pkin(7) + qJ(2);
t78 = sin(t82);
t79 = cos(t82);
t36 = t78 * t154 + t79 * t90;
t38 = -t79 * t154 + t78 * t90;
t170 = -g(1) * t38 + g(2) * t36;
t163 = g(3) * t89;
t169 = t86 * t163 + t170;
t130 = t86 * qJDD(2);
t167 = pkin(4) * t89;
t168 = t86 * pkin(4) * t112 + t130 * t167 + qJDD(5);
t165 = g(1) * t78;
t73 = g(2) * t79;
t129 = t87 * qJDD(2);
t65 = -qJDD(4) + t129;
t162 = t65 * pkin(4);
t144 = qJD(2) * t86;
t116 = qJ(5) * t144;
t43 = t56 * qJD(2) + qJD(3);
t132 = qJ(3) * qJD(2);
t54 = t86 * qJD(1) + t87 * t132;
t17 = t90 * t43 - t89 * t54;
t14 = -t90 * t116 + t17;
t134 = t87 * qJD(2);
t67 = -qJD(4) + t134;
t5 = -t67 * pkin(4) + t14;
t161 = -t14 + t5;
t77 = t87 * qJD(1);
t53 = t86 * t132 - t77;
t160 = t53 * t86;
t159 = t65 * t87;
t158 = t79 * t87;
t157 = t79 * t89;
t80 = t86 ^ 2;
t91 = qJD(2) ^ 2;
t156 = t80 * t91;
t155 = t86 * (-qJ(5) - pkin(6));
t137 = qJD(4) * t90;
t141 = qJD(3) * t87;
t152 = t56 * t137 + t90 * t141;
t151 = t79 * pkin(2) + t78 * qJ(3);
t81 = t87 ^ 2;
t150 = t80 + t81;
t83 = t89 ^ 2;
t84 = t90 ^ 2;
t149 = -t83 - t84;
t148 = t83 - t84;
t147 = qJ(3) * t89;
t146 = qJ(5) * t86;
t52 = (qJ(3) + t167) * t86;
t145 = qJD(2) * t52;
t143 = qJD(2) * t89;
t142 = qJD(2) * t90;
t18 = t89 * t43 + t90 * t54;
t140 = qJD(4) * t18;
t139 = qJD(4) * t54;
t138 = qJD(4) * t89;
t136 = qJD(5) * t86;
t135 = qJDD(2) * pkin(2);
t29 = qJD(5) - t77 + t145;
t133 = qJD(5) + t29;
t131 = qJDD(2) * t90;
t127 = qJ(5) * qJDD(2);
t124 = qJD(2) * qJD(5);
t123 = t89 * t156;
t122 = t87 * t147;
t121 = t90 * t146;
t120 = t86 * t143;
t119 = t67 * t138;
t118 = t53 * t144;
t117 = t89 * t141;
t115 = -t73 + t165;
t76 = t87 * qJDD(1);
t40 = t97 * t86 - t76;
t114 = -t40 * t87 - g(3);
t113 = t89 * t125;
t41 = t86 * qJDD(1) + t97 * t87;
t42 = t56 * qJDD(2) + qJDD(3);
t3 = t43 * t137 - t54 * t138 + t90 * t41 + t89 * t42;
t111 = t29 + t145;
t110 = t65 - t129;
t109 = t65 + t129;
t75 = qJDD(3) - t135;
t108 = t86 * t113;
t107 = t89 * t112;
t106 = -g(1) * t36 - g(2) * t38;
t37 = -t78 * t153 + t157;
t39 = t79 * t153 + t78 * t89;
t105 = -g(1) * t37 - g(2) * t39;
t104 = g(1) * t79 + g(2) * t78;
t15 = -t89 * t116 + t18;
t103 = t15 * t90 - t5 * t89;
t102 = -t15 * t89 - t5 * t90;
t101 = t17 * t89 - t18 * t90;
t100 = t40 * t86 + t41 * t87;
t99 = t54 * t87 + t160;
t98 = -t135 + t75 + t73;
t16 = t40 + t168;
t49 = (pkin(4) * t137 + qJD(3)) * t86;
t96 = qJD(2) * t49 + qJDD(2) * t52 + t16;
t32 = t90 * t42;
t94 = t32 + qJ(5) * t108 + (-qJD(4) * t43 - t41) * t89;
t93 = -t67 ^ 2 - t156;
t92 = g(3) * t86 * t90 + g(1) * t39 - g(2) * t37 - t3;
t4 = -t89 * t41 - t140 + t32;
t85 = qJDD(1) - g(3);
t74 = t90 * pkin(4) + pkin(3);
t71 = t79 * qJ(3);
t63 = t90 * t130;
t60 = t86 * t165;
t58 = t90 * t123;
t55 = t87 * t108;
t51 = t90 * t56;
t48 = t86 * t67 * t142;
t46 = t148 * t156;
t45 = t149 * t130;
t34 = (qJDD(2) * t84 - 0.2e1 * t107) * t80;
t33 = (qJDD(2) * t83 + 0.2e1 * t107) * t80;
t26 = t51 - t122;
t25 = -t89 * t146 + t27;
t24 = 0.2e1 * (t148 * t125 - t89 * t131) * t80;
t23 = -t48 - t171;
t22 = t63 + (-qJD(4) - t67) * t120;
t21 = -t117 - t172;
t20 = -qJD(4) * t122 + t152;
t19 = -t121 + t51 + (-pkin(4) - t147) * t87;
t13 = t89 * t65 + t93 * t90;
t12 = -t90 * t65 + t93 * t89;
t11 = -t117 - t90 * t136 + (-t66 + (-t56 + t146) * t89) * qJD(4);
t10 = -t89 * t136 + (-t121 - t122) * qJD(4) + t152;
t9 = (t110 * t89 + (t67 - t134) * t137) * t86;
t8 = (t109 * t89 + (t67 + t134) * t137) * t86;
t7 = t55 + (-t109 * t90 + t119) * t86;
t6 = t55 + (t110 * t90 - t119) * t86;
t2 = (-t95 * qJ(5) - t89 * t124) * t86 + t3;
t1 = -t162 + (-t139 + (-t124 - t127) * t86) * t90 + t94;
t28 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t86 + t114, 0, 0, 0, 0, 0, 0, t9, t6, 0, (t3 * t90 - t4 * t89 + (-t17 * t90 - t18 * t89) * qJD(4)) * t86 + t114, 0, 0, 0, 0, 0, 0, t9, t6, 0, -t16 * t87 - g(3) + (t102 * qJD(4) - t1 * t89 + t2 * t90) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t115, t104, 0, 0, t80 * qJDD(2), 0.2e1 * t86 * t129, 0, t81 * qJDD(2), 0, 0, (-t98 + t165) * t87, t98 * t86 - t60, t97 * t150 + t100 - t104, -t75 * pkin(2) - g(1) * (-t78 * pkin(2) + t71) - g(2) * t151 + t99 * qJD(3) + t100 * qJ(3), t34, t24, t7, t33, t8, t159, -t21 * t67 - t26 * t65 - t4 * t87 + (t53 * t137 + t40 * t89) * t86 + (t95 * qJ(3) + t89 * t126) * t80 + t105, t20 * t67 + t27 * t65 + t3 * t87 + (-t53 * t138 + t40 * t90) * t86 + (t90 * t126 + (-t113 + t131) * qJ(3)) * t80 + t106, t60 + (-t73 + (-t140 - qJDD(2) * t26 - t4 + (-t21 - t172) * qJD(2)) * t90 + (qJD(4) * t17 - qJDD(2) * t27 - t3 + (qJD(4) * t26 - t20) * qJD(2)) * t89) * t86, t3 * t27 + t18 * t20 + t4 * t26 + t17 * t21 - g(1) * t71 - g(2) * (pkin(3) * t158 + t151) + (-pkin(6) * t73 + t40 * qJ(3) + t53 * qJD(3)) * t86 - t56 * t165, t34, t24, t7, t33, t8, t159, -t1 * t87 - t11 * t67 - t19 * t65 + (t111 * t137 + t96 * t89) * t86 + t105, t10 * t67 + t2 * t87 + t25 * t65 + (-t111 * t138 + t96 * t90) * t86 + t106, t60 + (-t73 + (-qJD(4) * t15 - qJDD(2) * t19 - t1 + (-qJD(4) * t25 - t11) * qJD(2)) * t90 + (qJD(4) * t5 - qJDD(2) * t25 - t2 + (qJD(4) * t19 - t10) * qJD(2)) * t89) * t86, t2 * t25 + t15 * t10 + t1 * t19 + t5 * t11 + t16 * t52 + t29 * t49 - g(1) * (pkin(4) * t157 + t71) - g(2) * (-t79 * t155 + t74 * t158 + t151) + (-g(1) * (-t74 * t87 - pkin(2) + t155) - g(2) * t167) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t130, -t150 * t91, -t99 * qJD(2) - t115 + t75, 0, 0, 0, 0, 0, 0, t12, t13, t45, t3 * t89 + t4 * t90 - t101 * qJD(4) + (t101 * t87 - t160) * qJD(2) - t115, 0, 0, 0, 0, 0, 0, t12, t13, t45, t1 * t90 + t2 * t89 + t103 * qJD(4) + (-t103 * t87 - t29 * t86) * qJD(2) - t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t46, t22, -t58, t23, -t65, -t90 * t118 - t18 * t67 + t169 + t4, t89 * t118 - t17 * t67 + t92, 0, 0, t58, -t46, t22, -t58, t23, -t65, -0.2e1 * t162 - t15 * t67 + (-pkin(4) * t123 - t139 + (-t133 * qJD(2) - t127) * t86) * t90 + t94 + t169, -t84 * pkin(4) * t156 - t14 * t67 + (t89 * t127 + (qJ(5) * t137 + t133 * t89) * qJD(2)) * t86 + t92, (-pkin(4) * t131 + (pkin(4) * qJD(4) - t161) * t143) * t86, t161 * t15 + (t1 + (-t29 * t142 + t163) * t86 + t170) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 + t171, t63 + (-qJD(4) + t67) * t120, t149 * t156, g(3) * t87 - t76 + (t128 + (qJD(3) - t102) * qJD(2) - t104) * t86 + t168;];
tau_reg = t28;
