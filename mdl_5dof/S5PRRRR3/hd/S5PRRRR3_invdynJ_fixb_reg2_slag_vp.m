% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:24
% EndTime: 2019-12-05 17:06:26
% DurationCPUTime: 0.84s
% Computational Cost: add. (1593->186), mult. (2470->234), div. (0->0), fcn. (1355->12), ass. (0->122)
t70 = pkin(9) + qJ(2);
t68 = qJ(3) + t70;
t57 = qJ(4) + t68;
t52 = cos(t57);
t146 = g(2) * t52;
t69 = qJDD(2) + qJDD(3);
t63 = qJDD(4) + t69;
t144 = t63 * pkin(4);
t79 = cos(qJ(4));
t123 = qJD(4) * t79;
t77 = sin(qJ(3));
t115 = t77 * t123;
t121 = qJDD(2) * t77;
t76 = sin(qJ(4));
t124 = qJD(4) * t76;
t80 = cos(qJ(3));
t125 = qJD(3) * t80;
t127 = pkin(2) * qJD(2);
t117 = t77 * t127;
t142 = t80 * pkin(2);
t61 = qJDD(2) * t142;
t24 = t69 * pkin(3) - qJD(3) * t117 + t61;
t21 = t79 * t24;
t71 = qJD(2) + qJD(3);
t34 = t71 * pkin(3) + t80 * t127;
t8 = -pkin(2) * ((t76 * t125 + t115) * qJD(2) + t76 * t121) - t34 * t124 + t21;
t6 = -t144 - t8;
t114 = -t6 - t146;
t19 = t79 * t117 + t76 * t34;
t66 = qJD(4) + t71;
t17 = t66 * pkin(8) + t19;
t78 = cos(qJ(5));
t135 = t78 * t17;
t75 = sin(qJ(5));
t10 = t75 * qJD(1) + t135;
t9 = t78 * qJD(1) - t75 * t17;
t126 = t9 * qJD(5);
t106 = qJD(4) * t117;
t110 = qJD(2) * t125;
t136 = t77 * t79;
t53 = pkin(2) * t136;
t102 = -t79 * pkin(2) * t110 - qJDD(2) * t53 - t34 * t123 + (t106 - t24) * t76;
t5 = t63 * pkin(8) - t102;
t2 = t75 * qJDD(1) + t78 * t5 + t126;
t67 = t78 * qJDD(1);
t3 = -t10 * qJD(5) - t75 * t5 + t67;
t86 = t2 * t78 - t3 * t75 + (-t10 * t75 - t78 * t9) * qJD(5);
t51 = sin(t57);
t131 = -g(1) * t52 - g(2) * t51;
t137 = t76 * t77;
t95 = t79 * t80 - t137;
t28 = t95 * t127;
t151 = pkin(3) * t123 - t28;
t44 = g(1) * t51;
t150 = t44 - t146;
t55 = sin(t68);
t48 = g(1) * t55;
t56 = cos(t68);
t149 = -g(2) * t56 + t48;
t96 = t76 * t80 + t136;
t27 = t96 * t127;
t58 = t76 * pkin(3) + pkin(8);
t59 = -t79 * pkin(3) - pkin(4);
t81 = qJD(5) ^ 2;
t148 = -(pkin(3) * t124 - t27) * t66 - t58 * t81 - t59 * t63;
t143 = t66 * pkin(4);
t60 = pkin(3) + t142;
t11 = t60 * t123 + (t95 * qJD(3) - t77 * t124) * pkin(2);
t141 = t11 * t66;
t12 = t60 * t124 + (t96 * qJD(3) + t115) * pkin(2);
t140 = t12 * t66;
t18 = -t76 * t117 + t79 * t34;
t139 = t18 * t66;
t138 = t19 * t66;
t134 = t78 * t63;
t16 = -t18 - t143;
t133 = t16 * qJD(5) * t75 + t78 * t44;
t132 = t52 * pkin(4) + t51 * pkin(8);
t130 = g(1) * t56 + g(2) * t55;
t30 = t76 * t60 + t53;
t72 = t75 ^ 2;
t73 = t78 ^ 2;
t129 = t72 - t73;
t128 = t72 + t73;
t122 = qJD(5) * t78;
t74 = qJDD(1) - g(3);
t120 = -t114 * t75 + t16 * t122;
t62 = t66 ^ 2;
t119 = t75 * t62 * t78;
t50 = pkin(3) * t56;
t118 = t50 + t132;
t112 = -t51 * pkin(4) + t52 * pkin(8);
t111 = t128 * t63;
t109 = qJD(2) * (-qJD(3) + t71);
t108 = qJD(3) * (-qJD(2) - t71);
t107 = t75 * t66 * t122;
t105 = t61 + t149;
t103 = g(1) * t112;
t64 = sin(t70);
t101 = -pkin(2) * t64 - pkin(3) * t55;
t65 = cos(t70);
t100 = g(1) * t64 - g(2) * t65;
t98 = t10 * t78 - t9 * t75;
t29 = -pkin(2) * t137 + t79 * t60;
t94 = t102 - t131;
t93 = pkin(8) * t81 - t138 - t144;
t25 = -pkin(4) - t29;
t26 = pkin(8) + t30;
t91 = t25 * t63 + t26 * t81 + t140;
t90 = -pkin(8) * qJDD(5) + (t18 - t143) * qJD(5);
t89 = -qJDD(5) * t26 + (t25 * t66 - t11) * qJD(5);
t87 = -qJD(1) * qJD(5) - t16 * t66 - t131 - t5;
t85 = t131 + t86;
t84 = -qJDD(5) * t58 + (t59 * t66 - t151) * qJD(5);
t83 = t8 + t150;
t54 = pkin(2) * t65;
t36 = qJDD(5) * t78 - t81 * t75;
t35 = qJDD(5) * t75 + t81 * t78;
t23 = t73 * t63 - 0.2e1 * t107;
t22 = t72 * t63 + 0.2e1 * t107;
t15 = -0.2e1 * t129 * t66 * qJD(5) + 0.2e1 * t75 * t134;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, t36, -t35, 0, t98 * qJD(5) + t2 * t75 + t3 * t78 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t100, g(1) * t65 + g(2) * t64, 0, 0, 0, 0, 0, 0, 0, t69, (t77 * t108 + t69 * t80) * pkin(2) + t105, ((-qJDD(2) - t69) * t77 + t80 * t108) * pkin(2) + t130, 0, (t100 + (t77 ^ 2 + t80 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), 0, 0, 0, 0, 0, t63, t29 * t63 - t140 + t83, -t30 * t63 - t141 + t94, 0, -t102 * t30 + t19 * t11 + t8 * t29 - t18 * t12 - g(1) * t101 - g(2) * (t50 + t54), t22, t15, t35, t23, t36, 0, t89 * t75 + (t114 - t91) * t78 + t133, t89 * t78 + (t91 - t44) * t75 + t120, t26 * t111 + t128 * t141 + t85, t6 * t25 + t16 * t12 - g(1) * (t101 + t112) - g(2) * (t54 + t118) + t98 * t11 + t86 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t77 * pkin(2) * t109 + t105, (t80 * t109 - t121) * pkin(2) + t130, 0, 0, 0, 0, 0, 0, 0, t63, t27 * t66 + t21 + (pkin(3) * t63 - t106) * t79 + ((-pkin(3) * t66 - t34) * qJD(4) + (-t110 - t121) * pkin(2)) * t76 + t150, t28 * t66 + (-t66 * t123 - t63 * t76) * pkin(3) + t94, 0, t18 * t27 - t19 * t28 + (-t102 * t76 + t79 * t8 + (-t18 * t76 + t19 * t79) * qJD(4) + t149) * pkin(3), t22, t15, t35, t23, t36, 0, t84 * t75 + (t114 + t148) * t78 + t133, t84 * t78 + (-t148 - t44) * t75 + t120, t128 * t151 * t66 + t58 * t111 + t85, t6 * t59 - t16 * t27 - t103 - g(2) * t118 - t98 * t28 + (t48 + (t16 * t76 + t98 * t79) * qJD(4)) * pkin(3) + t86 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t83 + t138, t94 + t139, 0, 0, t22, t15, t35, t23, t36, 0, t90 * t75 + (t114 - t93) * t78 + t133, t90 * t78 + (t93 - t44) * t75 + t120, pkin(8) * t111 - t128 * t139 + t85, -t6 * pkin(4) + t86 * pkin(8) - g(2) * t132 - t16 * t19 - t98 * t18 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t129 * t62, t75 * t63, t119, t134, qJDD(5), -g(3) * t78 + t67 + (t10 - t135) * qJD(5) + t87 * t75, t126 + (qJD(5) * t17 - t74) * t75 + t87 * t78, 0, 0;];
tau_reg = t1;
