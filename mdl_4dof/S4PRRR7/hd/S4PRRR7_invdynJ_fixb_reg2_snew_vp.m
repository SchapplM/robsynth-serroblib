% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:51
% DurationCPUTime: 0.71s
% Computational Cost: add. (1604->165), mult. (3156->234), div. (0->0), fcn. (2252->10), ass. (0->111)
t117 = -g(3) + qJDD(1);
t84 = sin(pkin(4));
t85 = cos(pkin(4));
t119 = sin(pkin(8));
t120 = cos(pkin(8));
t99 = t119 * g(1) - t120 * g(2);
t98 = t85 * t99;
t132 = t84 * t117 + t98;
t95 = qJD(2) ^ 2;
t89 = sin(qJ(3));
t118 = qJD(2) * t89;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t53 = -t91 * qJD(3) + t88 * t118;
t55 = t88 * qJD(3) + t91 * t118;
t129 = t55 * t53;
t113 = qJD(2) * qJD(3);
t74 = t89 * t113;
t92 = cos(qJ(3));
t76 = t92 * qJDD(2);
t59 = t76 - t74;
t52 = -qJDD(4) + t59;
t100 = -t52 - t129;
t131 = t100 * t88;
t130 = t100 * t91;
t111 = t92 * t113;
t114 = t89 * qJDD(2);
t58 = t111 + t114;
t110 = -t91 * qJDD(3) + t88 * t58;
t71 = t92 * qJD(2) - qJD(4);
t24 = (qJD(4) + t71) * t55 + t110;
t50 = t53 ^ 2;
t51 = t55 ^ 2;
t69 = t71 ^ 2;
t128 = t71 * t88;
t127 = t71 * t91;
t106 = -t92 * pkin(3) - t89 * pkin(7);
t63 = -t120 * g(1) - t119 * g(2);
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t37 = t132 * t90 + t93 * t63;
t31 = -t95 * pkin(2) + qJDD(2) * pkin(6) + t37;
t108 = t95 * t106 + t31;
t97 = t85 * t117 - t84 * t99;
t45 = t92 * t97;
t94 = qJD(3) ^ 2;
t13 = -qJDD(3) * pkin(3) - t94 * pkin(7) + t108 * t89 - t45;
t126 = t88 * t13;
t33 = t52 - t129;
t125 = t88 * t33;
t70 = t89 * t95 * t92;
t64 = qJDD(3) + t70;
t124 = t89 * t64;
t123 = t91 * t13;
t122 = t91 * t33;
t65 = qJDD(3) - t70;
t121 = t92 * t65;
t116 = qJD(4) - t71;
t112 = t92 * t129;
t96 = t89 * t97;
t14 = -t94 * pkin(3) + qJDD(3) * pkin(7) + t108 * t92 + t96;
t103 = -t59 + t74;
t104 = t58 + t111;
t107 = -t132 * t93 + t90 * t63;
t30 = -qJDD(2) * pkin(2) - t95 * pkin(6) + t107;
t16 = t103 * pkin(3) - t104 * pkin(7) + t30;
t5 = t88 * t14 - t91 * t16;
t6 = t91 * t14 + t88 * t16;
t3 = t88 * t5 + t91 * t6;
t22 = t89 * t31 - t45;
t23 = t92 * t31 + t96;
t8 = t89 * t22 + t92 * t23;
t105 = t91 * t5 - t88 * t6;
t102 = -t88 * qJDD(3) - t91 * t58;
t101 = -pkin(2) + t106;
t39 = -t53 * qJD(4) - t102;
t81 = t92 ^ 2;
t80 = t89 ^ 2;
t79 = t81 * t95;
t77 = t80 * t95;
t68 = -t79 - t94;
t67 = -t77 - t94;
t62 = t77 + t79;
t61 = (t80 + t81) * qJDD(2);
t60 = t76 - 0.2e1 * t74;
t57 = 0.2e1 * t111 + t114;
t48 = t53 * t71;
t47 = -t51 + t69;
t46 = t50 - t69;
t44 = -t89 * t67 - t121;
t43 = t92 * t68 - t124;
t42 = t51 - t50;
t41 = -t51 - t69;
t40 = -t69 - t50;
t38 = -t55 * qJD(4) - t110;
t32 = t50 + t51;
t29 = t116 * t53 + t102;
t28 = t39 - t48;
t27 = t39 + t48;
t25 = -t116 * t55 - t110;
t21 = -t88 * t41 + t122;
t20 = t91 * t41 + t125;
t18 = t91 * t40 - t131;
t17 = t88 * t40 + t130;
t12 = -t24 * t91 + t88 * t28;
t11 = -t24 * t88 - t91 * t28;
t10 = t92 * t21 - t89 * t29;
t9 = t92 * t18 - t89 * t25;
t7 = t92 * t12 - t89 * t32;
t1 = t89 * t13 + t92 * t3;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, (qJDD(2) * t93 - t90 * t95) * t84, (-qJDD(2) * t90 - t93 * t95) * t84, 0, t85 ^ 2 * t117 + (-t107 * t93 + t90 * t37 - t98) * t84, 0, 0, 0, 0, 0, 0, t85 * (t92 * t64 + t89 * t68) + (t90 * t43 + t93 * t60) * t84, t85 * (-t89 * t65 + t92 * t67) + (t90 * t44 - t93 * t57) * t84, (t61 * t90 + t62 * t93) * t84, t85 * (-t92 * t22 + t89 * t23) + (-t93 * t30 + t90 * t8) * t84, 0, 0, 0, 0, 0, 0, t85 * (t89 * t18 + t92 * t25) + (-t93 * t17 + t90 * t9) * t84, t85 * (t89 * t21 + t92 * t29) + (t90 * t10 - t93 * t20) * t84, t85 * (t89 * t12 + t92 * t32) + (-t93 * t11 + t90 * t7) * t84, t85 * (-t92 * t13 + t89 * t3) + (t90 * t1 + t105 * t93) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t107, -t37, 0, 0, t104 * t89, t92 * t57 + t89 * t60, t124 + t92 * (-t77 + t94), -t103 * t92, t89 * (t79 - t94) + t121, 0, pkin(2) * t60 + pkin(6) * t43 - t92 * t30, -pkin(2) * t57 + pkin(6) * t44 + t89 * t30, pkin(2) * t62 + pkin(6) * t61 + t8, -pkin(2) * t30 + pkin(6) * t8, t89 * (t55 * t128 + t91 * t39) - t112, t89 * (t91 * t25 - t88 * t27) - t92 * t42, t89 * (-t88 * t47 + t130) - t92 * t28, t89 * (-t53 * t127 - t88 * t38) + t112, t89 * (t91 * t46 + t125) + t92 * t24, t92 * t52 + t89 * (t53 * t91 - t55 * t88) * t71, t89 * (-pkin(7) * t17 + t126) + t92 * (-pkin(3) * t17 + t5) - pkin(2) * t17 + pkin(6) * t9, t89 * (-pkin(7) * t20 + t123) + t92 * (-pkin(3) * t20 + t6) - pkin(2) * t20 + pkin(6) * t10, pkin(6) * t7 + t101 * t11 + t89 * t105, pkin(6) * t1 - t101 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t77 - t79, t114, t70, t76, qJDD(3), -t22, -t23, 0, 0, -t55 * t127 + t88 * t39, t88 * t25 + t91 * t27, t91 * t47 + t131, -t53 * t128 + t91 * t38, t88 * t46 - t122, (t53 * t88 + t55 * t91) * t71, pkin(3) * t25 + pkin(7) * t18 - t123, pkin(3) * t29 + pkin(7) * t21 + t126, pkin(3) * t32 + pkin(7) * t12 + t3, -pkin(3) * t13 + pkin(7) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t42, t28, -t129, -t24, -t52, -t5, -t6, 0, 0;];
tauJ_reg = t2;
