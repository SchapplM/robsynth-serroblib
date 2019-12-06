% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:20
% EndTime: 2019-12-05 15:07:21
% DurationCPUTime: 0.92s
% Computational Cost: add. (907->196), mult. (2065->244), div. (0->0), fcn. (1540->10), ass. (0->126)
t123 = qJD(1) * qJD(3);
t80 = sin(qJ(3));
t116 = t80 * t123;
t82 = cos(qJ(3));
t126 = qJDD(1) * t82;
t160 = qJDD(1) * t80 + t82 * t123;
t74 = sin(pkin(8));
t76 = cos(pkin(8));
t119 = -t74 * t126 - t160 * t76;
t95 = -t116 * t74 - t119;
t10 = qJDD(3) * pkin(6) + t95;
t161 = qJD(2) * qJD(4) + t10;
t71 = pkin(8) + qJ(3);
t64 = cos(t71);
t154 = g(3) * t64;
t75 = sin(pkin(7));
t77 = cos(pkin(7));
t105 = g(1) * t77 + g(2) * t75;
t63 = sin(t71);
t98 = t105 * t63;
t89 = -t154 + t98;
t38 = t74 * t82 + t76 * t80;
t32 = t38 * qJD(1);
t24 = qJD(3) * pkin(6) + t32;
t107 = qJ(5) * qJD(3) + t24;
t81 = cos(qJ(4));
t101 = t107 * t81;
t118 = -g(1) * t75 + g(2) * t77;
t155 = g(3) * t63;
t91 = t105 * t64 + t155;
t79 = sin(qJ(4));
t72 = t79 ^ 2;
t159 = pkin(4) * t72;
t153 = g(3) * t79;
t152 = t81 * pkin(4);
t151 = t75 * t79;
t150 = t75 * t81;
t149 = t77 * t79;
t148 = t77 * t81;
t147 = t80 * t74;
t144 = t82 * t76;
t57 = qJD(1) * t147;
t31 = qJD(1) * t144 - t57;
t146 = t81 * t31;
t84 = qJD(3) ^ 2;
t145 = t81 * t84;
t143 = qJ(5) + pkin(6);
t69 = t81 * qJD(2);
t13 = -t107 * t79 + t69;
t137 = qJD(4) * pkin(4);
t12 = t13 + t137;
t142 = -t13 + t12;
t141 = qJD(4) * t146 + t153 * t64;
t73 = t81 ^ 2;
t140 = t72 - t73;
t139 = t72 + t73;
t138 = qJD(3) * pkin(3);
t37 = -t144 + t147;
t136 = qJD(3) * t37;
t135 = qJD(3) * t81;
t134 = qJD(4) * t79;
t133 = qJDD(3) * pkin(3);
t132 = qJDD(4) * pkin(4);
t131 = t32 * qJD(3);
t130 = t136 * qJD(3);
t129 = t79 * qJD(2);
t62 = pkin(3) + t152;
t17 = -qJD(3) * t62 + qJD(5) - t31;
t128 = -qJD(5) - t17;
t125 = qJDD(3) * t62;
t65 = t79 * qJDD(3);
t67 = t81 * qJDD(3);
t124 = qJ(5) * qJDD(3);
t121 = qJD(3) * qJD(4);
t120 = qJD(3) * qJD(5);
t117 = t139 * t31;
t114 = t79 * t121;
t113 = t81 * t121;
t6 = t79 * qJDD(2) - t24 * t134 + t161 * t81;
t112 = qJD(4) * t143;
t111 = t81 * t131 + t31 * t134 + (g(1) * t148 + g(2) * t150) * t63;
t110 = (t116 - t126) * t76 + t160 * t74;
t8 = pkin(4) * t114 + qJDD(5) + t110 - t125;
t109 = t8 - t125;
t108 = t139 * qJDD(3);
t106 = t79 * t113;
t14 = t129 + t101;
t104 = t12 * t79 - t14 * t81;
t15 = -t24 * t79 + t69;
t16 = t24 * t81 + t129;
t103 = t15 * t79 - t16 * t81;
t34 = t38 * qJD(3);
t100 = -qJD(3) * t34 - qJDD(3) * t37;
t11 = t110 - t133;
t83 = qJD(4) ^ 2;
t99 = pkin(6) * t83 + t11 - t133;
t68 = t81 * qJDD(2);
t97 = -g(1) * (-t149 * t64 + t150) - g(2) * (-t151 * t64 - t148) + t63 * t153 + t68;
t96 = -t124 - t161;
t94 = t38 * t83 - t100;
t23 = -t31 - t138;
t93 = -pkin(6) * qJDD(4) + (t23 - t138) * qJD(4);
t90 = 0.2e1 * qJD(4) * t136 - qJDD(4) * t38;
t88 = -g(1) * (-t148 * t64 - t151) - g(2) * (-t150 * t64 + t149) + t81 * t155 - t6;
t87 = -t98 - t131;
t7 = -t16 * qJD(4) - t79 * t10 + t68;
t86 = t6 * t81 - t7 * t79 + (-t15 * t81 - t16 * t79) * qJD(4);
t85 = t86 - t91;
t58 = t79 * t145;
t46 = t143 * t81;
t45 = t143 * t79;
t44 = t140 * t84;
t43 = qJDD(4) * t81 - t79 * t83;
t42 = qJDD(4) * t79 + t81 * t83;
t41 = qJDD(2) + t118;
t36 = qJDD(3) * t73 - 0.2e1 * t106;
t35 = qJDD(3) * t72 + 0.2e1 * t106;
t30 = -t79 * qJD(5) - t112 * t81;
t29 = t81 * qJD(5) - t112 * t79;
t19 = -0.2e1 * t121 * t140 + 0.2e1 * t67 * t79;
t5 = t108 * t38 - t130 * t139;
t4 = t81 * t120 + (-t114 + t67) * qJ(5) + t6;
t3 = t132 + t68 - qJD(4) * t101 + (t96 - t120) * t79;
t2 = t79 * t90 - t81 * t94;
t1 = t79 * t94 + t81 * t90;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t74 ^ 2 + t76 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t100, -qJDD(3) * t38 + t130, 0, t110 * t37 - t136 * t32 - t31 * t34 + t38 * t95 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t5, t103 * t136 + t11 * t37 + t23 * t34 + t38 * t86 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t5, t17 * t34 + t8 * t37 - g(3) + t104 * t136 + (-t3 * t79 + t4 * t81 + (-t12 * t81 - t14 * t79) * qJD(4)) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t43, -t42, 0, -qJD(4) * t103 + t6 * t79 + t7 * t81 + t118, 0, 0, 0, 0, 0, 0, t43, -t42, 0, -qJD(4) * t104 + t3 * t81 + t4 * t79 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t110 + t89 + t131, (t31 + t57) * qJD(3) + t91 + t119, 0, 0, t35, t19, t42, t36, t43, 0, t93 * t79 + (-t99 - t154) * t81 + t111, t93 * t81 + (t87 + t99) * t79 + t141, pkin(6) * t108 - qJD(3) * t117 + t85, -t23 * t32 + t103 * t31 + (-t11 + t89) * pkin(3) + t85 * pkin(6), t35, t19, t42, t36, t43, 0, -t45 * qJDD(4) + (-t109 - t154) * t81 + (t30 + (t17 + (-t62 - t152) * qJD(3)) * t79) * qJD(4) + t111, -t46 * qJDD(4) + (t17 * t81 - t29 + (-t62 * t81 + t159) * qJD(3)) * qJD(4) + (t109 + t87) * t79 + t141, (-qJD(4) * t12 + qJDD(3) * t46 + t4) * t81 + (-t14 * qJD(4) + qJDD(3) * t45 - t3) * t79 + (t29 * t81 - t30 * t79 - t117 + (t45 * t81 - t46 * t79) * qJD(4)) * qJD(3) - t91, t4 * t46 - t3 * t45 - t8 * t62 - g(3) * (t143 * t63 + t62 * t64) + (pkin(4) * t134 - t32) * t17 + (t29 - t146) * t14 + (t31 * t79 + t30) * t12 + t105 * (-t143 * t64 + t62 * t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t44, t65, t58, t67, qJDD(4), (-qJD(3) * t23 - t10) * t79 + t97, qJD(4) * t15 - t135 * t23 + t88, 0, 0, -t58, t44, t65, t58, t67, qJDD(4), 0.2e1 * t132 + (t14 - t101) * qJD(4) + (pkin(4) * t145 + qJD(3) * t128 + t96) * t79 + t97, -t84 * t159 - t81 * t124 + t13 * qJD(4) + (qJ(5) * t134 + t128 * t81) * qJD(3) + t88, -pkin(4) * t65 + (-t137 + t142) * t135, t142 * t14 + (t3 + t118 * t81 + (-qJD(3) * t17 + t91) * t79) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 + 0.2e1 * t114, t65 + 0.2e1 * t113, -t139 * t84, qJD(3) * t104 + t8 - t89;];
tau_reg = t9;
