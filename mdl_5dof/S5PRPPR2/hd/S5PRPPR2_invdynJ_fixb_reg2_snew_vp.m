% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:42
% EndTime: 2019-12-05 15:24:47
% DurationCPUTime: 0.89s
% Computational Cost: add. (2075->159), mult. (4080->252), div. (0->0), fcn. (2812->10), ass. (0->101)
t121 = qJD(2) ^ 2;
t113 = sin(pkin(7));
t114 = cos(pkin(7));
t100 = -t114 * g(1) - t113 * g(2);
t109 = -g(3) + qJDD(1);
t90 = sin(qJ(2));
t92 = cos(qJ(2));
t53 = t92 * t100 + t90 * t109;
t51 = -t121 * pkin(2) + t53;
t86 = sin(pkin(8));
t88 = cos(pkin(8));
t52 = -t90 * t100 + t92 * t109;
t98 = qJDD(2) * pkin(2) + t52;
t33 = t88 * t51 + t86 * t98;
t126 = -t121 * pkin(3) + qJDD(2) * qJ(4) + (2 * qJD(2) * qJD(4)) + t33;
t85 = sin(pkin(9));
t87 = cos(pkin(9));
t89 = sin(qJ(5));
t91 = cos(qJ(5));
t62 = (t85 * t89 - t87 * t91) * qJD(2);
t101 = t85 * t91 + t87 * t89;
t64 = t101 * qJD(2);
t44 = t64 * t62;
t122 = qJDD(5) - t44;
t128 = t122 * t89;
t127 = t122 * t91;
t94 = t85 ^ 2;
t96 = t87 ^ 2;
t125 = t94 + t96;
t124 = t121 * t87;
t69 = -t113 * g(1) + t114 * g(2) + qJDD(3);
t66 = t87 * t69;
t123 = t66 + (pkin(4) * t124 - pkin(6) * qJDD(2) - t126) * t85;
t73 = t125 * t121;
t59 = t62 ^ 2;
t60 = t64 ^ 2;
t106 = t87 * qJDD(2);
t110 = t96 * t121;
t21 = t126 * t87 + t85 * t69;
t18 = -pkin(4) * t110 + pkin(6) * t106 + t21;
t8 = -t91 * t123 + t89 * t18;
t9 = t123 * t89 + t91 * t18;
t3 = -t91 * t8 + t89 * t9;
t120 = t85 * t3;
t32 = -t86 * t51 + t88 * t98;
t84 = qJDD(2) * pkin(3);
t28 = -t121 * qJ(4) + qJDD(4) - t32 - t84;
t22 = -pkin(4) * t106 + t28 + (-t121 * t94 - t110) * pkin(6);
t119 = t89 * t22;
t38 = qJDD(5) + t44;
t118 = t89 * t38;
t117 = t91 * t22;
t116 = t91 * t38;
t112 = t62 * qJD(5);
t111 = t64 * qJD(5);
t108 = t85 * qJDD(2);
t107 = t86 * qJDD(2);
t105 = t88 * qJDD(2);
t4 = t89 * t8 + t91 * t9;
t20 = t126 * t85 - t66;
t11 = t85 * t20 + t87 * t21;
t102 = -t28 + t84;
t34 = t91 * t106 - t89 * t108;
t61 = t101 * qJDD(2);
t93 = qJD(5) ^ 2;
t82 = t96 * qJDD(2);
t81 = t94 * qJDD(2);
t72 = -t88 * t121 - t107;
t71 = -t86 * t121 + t105;
t70 = t82 + t81;
t68 = t125 * t124;
t67 = t85 * t73;
t56 = -t60 - t93;
t55 = -t60 + t93;
t54 = t59 - t93;
t50 = t87 * t105 - t86 * t68;
t49 = -t85 * t105 + t86 * t67;
t45 = t86 * t70 + t88 * t73;
t43 = t61 - t112;
t42 = t61 - 0.2e1 * t112;
t41 = t34 - t111;
t40 = -t34 + 0.2e1 * t111;
t36 = -t93 - t59;
t35 = -t59 - t60;
t30 = -t89 * t56 - t116;
t29 = t91 * t56 - t118;
t26 = t91 * t34 + t89 * t61;
t25 = t89 * t34 - t91 * t61;
t24 = t91 * t36 - t128;
t23 = t89 * t36 + t127;
t16 = t88 * t32 + t86 * t33;
t15 = -t85 * t29 + t87 * t30;
t14 = -t85 * t25 + t87 * t26;
t13 = -t85 * t23 + t87 * t24;
t12 = t86 * t15 - t88 * t42;
t10 = t86 * t13 - t88 * t40;
t6 = t86 * t14 - t88 * t35;
t5 = t86 * t11 - t88 * t28;
t2 = t87 * t4 - t120;
t1 = t86 * t2 - t88 * t22;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, t92 * qJDD(2) - t90 * t121, -t90 * qJDD(2) - t92 * t121, 0, t92 * t52 + t90 * t53, 0, 0, 0, 0, 0, 0, t92 * t71 + t90 * t72, -t90 * t71 + t92 * t72, 0, t90 * (-t86 * t32 + t88 * t33) + t92 * t16, 0, 0, 0, 0, 0, 0, t90 * (-t86 * t106 - t88 * t68) + t92 * t50, t90 * (t85 * t107 + t88 * t67) + t92 * t49, t90 * (t88 * t70 - t86 * t73) + t92 * t45, t90 * (t88 * t11 + t86 * t28) + t92 * t5, 0, 0, 0, 0, 0, 0, t90 * (t88 * t13 + t86 * t40) + t92 * t10, t90 * (t88 * t15 + t86 * t42) + t92 * t12, t90 * (t88 * t14 + t86 * t35) + t92 * t6, t90 * (t88 * t2 + t86 * t22) + t92 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t52, -t53, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t71 + t32, pkin(2) * t72 - t33, 0, pkin(2) * t16, t81, 0.2e1 * t85 * t106, 0, t82, 0, 0, pkin(2) * t50 - qJ(4) * t68 + t102 * t87, pkin(2) * t49 + qJ(4) * t67 - t102 * t85, pkin(2) * t45 + pkin(3) * t73 + qJ(4) * t70 + t11, pkin(2) * t5 - pkin(3) * t28 + qJ(4) * t11, t85 * (-t89 * t111 + t91 * t43) + t87 * (t91 * t111 + t89 * t43), t85 * (-t91 * t40 - t89 * t42) + t87 * (-t89 * t40 + t91 * t42), t85 * (-t89 * t55 + t127) + t87 * (t91 * t55 + t128), t85 * (t91 * t112 - t89 * t41) + t87 * (t89 * t112 + t91 * t41), t85 * (t91 * t54 - t118) + t87 * (t89 * t54 + t116), (t85 * (-t62 * t91 + t64 * t89) + t87 * (-t62 * t89 - t64 * t91)) * qJD(5), t85 * (-pkin(6) * t23 + t119) + t87 * (-pkin(4) * t40 + pkin(6) * t24 - t117) - pkin(3) * t40 + qJ(4) * t13 + pkin(2) * t10, t85 * (-pkin(6) * t29 + t117) + t87 * (-pkin(4) * t42 + pkin(6) * t30 + t119) - pkin(3) * t42 + qJ(4) * t15 + pkin(2) * t12, t85 * (-pkin(6) * t25 - t3) + t87 * (-pkin(4) * t35 + pkin(6) * t26 + t4) - pkin(3) * t35 + qJ(4) * t14 + pkin(2) * t6, -pkin(6) * t120 + t87 * (-pkin(4) * t22 + pkin(6) * t4) - pkin(3) * t22 + qJ(4) * t2 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t20 + t85 * t21, 0, 0, 0, 0, 0, 0, t87 * t23 + t85 * t24, t87 * t29 + t85 * t30, t87 * t25 + t85 * t26, t87 * t3 + t85 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, t108, -t73, t28, 0, 0, 0, 0, 0, 0, t40, t42, t35, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t60 - t59, t61, -t44, t34, qJDD(5), -t8, -t9, 0, 0;];
tauJ_reg = t7;
