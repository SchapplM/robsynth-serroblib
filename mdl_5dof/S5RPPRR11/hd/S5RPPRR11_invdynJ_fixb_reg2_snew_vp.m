% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:56
% DurationCPUTime: 0.65s
% Computational Cost: add. (1617->158), mult. (3059->187), div. (0->0), fcn. (1632->6), ass. (0->110)
t83 = cos(qJ(4));
t111 = qJD(1) * t83;
t79 = sin(qJ(5));
t82 = cos(qJ(5));
t55 = -t82 * qJD(4) + t79 * t111;
t57 = t79 * qJD(4) + t82 * t111;
t41 = t57 * t55;
t80 = sin(qJ(4));
t107 = t80 * qJDD(1);
t105 = qJD(1) * qJD(4);
t99 = t83 * t105;
t59 = -t99 - t107;
t53 = qJDD(5) - t59;
t124 = -t41 + t53;
t126 = t124 * t79;
t125 = t124 * t82;
t66 = t80 * qJD(1) + qJD(5);
t100 = t80 * t105;
t106 = t83 * qJDD(1);
t60 = -t100 + t106;
t98 = -t82 * qJDD(4) + t79 * t60;
t21 = (qJD(5) - t66) * t57 + t98;
t51 = t55 ^ 2;
t52 = t57 ^ 2;
t65 = t66 ^ 2;
t86 = qJD(1) ^ 2;
t123 = 2 * qJD(3);
t122 = t66 * t79;
t121 = t66 * t82;
t76 = t80 ^ 2;
t120 = t76 * t86;
t77 = t83 ^ 2;
t119 = t77 * t86;
t29 = t41 + t53;
t118 = t79 * t29;
t103 = t80 * t86 * t83;
t117 = t80 * (qJDD(4) + t103);
t116 = t82 * t29;
t114 = pkin(1) + qJ(3);
t75 = qJDD(1) * qJ(2);
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t95 = t84 * g(1) + t81 * g(2);
t92 = 0.2e1 * qJD(2) * qJD(1) - t95;
t89 = qJDD(3) + t92;
t45 = -t114 * t86 + t75 + t89;
t40 = -qJDD(1) * pkin(6) + t45;
t33 = t80 * g(3) + t83 * t40;
t85 = qJD(4) ^ 2;
t96 = pkin(4) * t80 - pkin(7) * t83;
t91 = t86 * t96;
t19 = qJDD(4) * pkin(4) + t85 * pkin(7) - t83 * t91 + t33;
t115 = t83 * t19;
t113 = qJ(2) - pkin(6);
t112 = t76 + t77;
t110 = qJDD(1) * pkin(1);
t108 = qJD(5) + t66;
t104 = t80 * t41;
t102 = qJD(1) * t123;
t78 = t86 * pkin(6);
t101 = t81 * g(1) - t84 * g(2);
t94 = -qJDD(2) + t101;
t88 = t86 * qJ(2) + t94;
t97 = t114 * qJDD(1);
t87 = t97 + t88;
t18 = -t59 * pkin(4) - t60 * pkin(7) - t78 + (t123 + (pkin(4) * t83 + pkin(7) * t80) * qJD(4)) * qJD(1) + t87;
t34 = t83 * g(3) - t80 * t40;
t20 = -t85 * pkin(4) + qJDD(4) * pkin(7) - t80 * t91 - t34;
t8 = -t82 * t18 + t79 * t20;
t9 = t79 * t18 + t82 * t20;
t3 = t79 * t8 + t82 * t9;
t2 = t79 * t9 - t82 * t8;
t14 = t83 * t33 - t80 * t34;
t93 = -t79 * qJDD(4) - t82 * t60;
t90 = t96 + t114;
t32 = -t55 * qJD(5) - t93;
t44 = t87 + t102;
t71 = 0.2e1 * t75;
t63 = t112 * t86;
t62 = t112 * qJDD(1);
t61 = -0.2e1 * t100 + t106;
t58 = 0.2e1 * t99 + t107;
t54 = t83 * (qJDD(4) - t103);
t49 = t88 + t110;
t48 = t66 * t55;
t47 = -t52 + t65;
t46 = t51 - t65;
t43 = -t117 + t83 * (-t85 - t119);
t42 = t80 * (-t85 - t120) + t54;
t39 = -t78 + t44;
t38 = t52 - t51;
t37 = -t52 - t65;
t35 = -t65 - t51;
t31 = -t57 * qJD(5) - t98;
t27 = t51 + t52;
t26 = t108 * t55 + t93;
t25 = t32 + t48;
t24 = t32 - t48;
t22 = -t108 * t57 - t98;
t16 = -t79 * t37 - t116;
t15 = t82 * t37 - t118;
t13 = t82 * t35 - t126;
t12 = t79 * t35 + t125;
t11 = -t21 * t82 + t79 * t25;
t10 = -t21 * t79 - t82 * t25;
t6 = t80 * t16 + t83 * t26;
t5 = t80 * t13 + t83 * t22;
t4 = t80 * t11 + t83 * t27;
t1 = t80 * t3 + t115;
t7 = [0, 0, 0, 0, 0, qJDD(1), t101, t95, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t94 - 0.2e1 * t110, t71 + t92, qJ(2) * (-t86 * pkin(1) + t75 + t92) + pkin(1) * t49, qJDD(1), 0, 0, 0, 0, 0, 0, t71 + t89, t102 + t94 + 0.2e1 * t97, qJ(2) * t45 + t114 * t44, (t60 - t100) * t83, -t83 * t58 - t80 * t61, t54 - t80 * (t85 - t119), (-t59 + t99) * t80, t83 * (-t85 + t120) - t117, 0, t113 * t42 + t114 * t58 + t80 * t39, t113 * t43 + t114 * t61 + t83 * t39, -t113 * t62 - t114 * t63 - t14, t113 * t14 + t114 * t39, t83 * (-t57 * t122 + t82 * t32) + t104, t83 * (t82 * t22 - t79 * t24) + t80 * t38, t83 * (-t79 * t47 + t125) + t80 * t25, t83 * (t55 * t121 - t79 * t31) - t104, t83 * (t82 * t46 - t118) - t80 * t21, t80 * t53 + t83 * (-t55 * t82 + t57 * t79) * t66, t113 * t5 - t79 * t115 + t90 * t12 - t80 * t8, t113 * t6 - t82 * t115 + t90 * t15 - t80 * t9, t90 * t10 + t113 * t4 - t83 * t2, t113 * t1 + t90 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t86, -t49, 0, 0, 0, 0, 0, 0, 0, -t86, -qJDD(1), -t44, 0, 0, 0, 0, 0, 0, -t58, -t61, t63, -t39, 0, 0, 0, 0, 0, 0, -t12, -t15, -t10, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t86, t45, 0, 0, 0, 0, 0, 0, t42, t43, -t62, t14, 0, 0, 0, 0, 0, 0, t5, t6, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, (-t76 + t77) * t86, t106, -t103, -t107, qJDD(4), t33, t34, 0, 0, t57 * t121 + t79 * t32, t79 * t22 + t82 * t24, t82 * t47 + t126, t55 * t122 + t82 * t31, t79 * t46 + t116, (-t55 * t79 - t57 * t82) * t66, pkin(4) * t22 + pkin(7) * t13 + t82 * t19, pkin(4) * t26 + pkin(7) * t16 - t79 * t19, pkin(4) * t27 + pkin(7) * t11 + t3, pkin(4) * t19 + pkin(7) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t38, t25, -t41, -t21, t53, -t8, -t9, 0, 0;];
tauJ_reg = t7;
