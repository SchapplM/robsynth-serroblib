% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:53
% DurationCPUTime: 0.87s
% Computational Cost: add. (2183->178), mult. (4141->263), div. (0->0), fcn. (2986->10), ass. (0->105)
t89 = sin(qJ(4));
t109 = qJD(3) * t89;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t92 = cos(qJ(4));
t53 = -t91 * t92 * qJD(3) + t88 * t109;
t55 = (t92 * t88 + t89 * t91) * qJD(3);
t41 = t55 * t53;
t80 = qJDD(4) + qJDD(5);
t124 = -t41 + t80;
t126 = t124 * t88;
t125 = t124 * t91;
t105 = qJD(3) * qJD(4);
t104 = t89 * t105;
t76 = t92 * qJDD(3);
t100 = t76 - t104;
t103 = t92 * t105;
t107 = t89 * qJDD(3);
t60 = t103 + t107;
t34 = -t53 * qJD(5) + t88 * t100 + t91 * t60;
t81 = qJD(4) + qJD(5);
t50 = t81 * t53;
t123 = -t50 + t34;
t108 = -g(3) + qJDD(1);
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t110 = sin(pkin(8));
t111 = cos(pkin(8));
t98 = -t111 * g(1) - t110 * g(2);
t46 = t85 * t108 + t86 * t98;
t90 = sin(qJ(3));
t93 = cos(qJ(3));
t97 = t86 * t108 - t85 * t98;
t32 = t93 * t46 + t90 * t97;
t95 = qJD(3) ^ 2;
t30 = -t95 * pkin(3) + qJDD(3) * pkin(6) + t32;
t62 = -t110 * g(1) + t111 * g(2) + qJDD(2);
t19 = t89 * t30 - t92 * t62;
t122 = -t19 + (-t60 + t103) * pkin(7);
t51 = t53 ^ 2;
t52 = t55 ^ 2;
t79 = t81 ^ 2;
t20 = t92 * t30 + t89 * t62;
t69 = qJD(4) * pkin(4) - pkin(7) * t109;
t83 = t92 ^ 2;
t78 = t83 * t95;
t14 = -pkin(4) * t78 + t100 * pkin(7) - qJD(4) * t69 + t20;
t72 = t89 * t95 * t92;
t106 = qJDD(4) + t72;
t96 = t106 * pkin(4) + t122;
t6 = t88 * t14 - t91 * t96;
t115 = t91 * t14;
t7 = t88 * t96 + t115;
t2 = -t91 * t6 + t88 * t7;
t121 = t89 * t2;
t120 = t81 * t88;
t119 = t81 * t91;
t31 = -t90 * t46 + t93 * t97;
t29 = -qJDD(3) * pkin(3) - t95 * pkin(6) - t31;
t15 = -t100 * pkin(4) - pkin(7) * t78 + t69 * t109 + t29;
t118 = t88 * t15;
t38 = t41 + t80;
t117 = t88 * t38;
t116 = t89 * t106;
t114 = t91 * t15;
t113 = t91 * t38;
t68 = qJDD(4) - t72;
t112 = t92 * t68;
t3 = t88 * t6 + t91 * t7;
t9 = t89 * t19 + t92 * t20;
t102 = -t91 * t100 + t88 * t60;
t61 = t76 - 0.2e1 * t104;
t99 = (-qJD(5) + t81) * t55 - t102;
t94 = qJD(4) ^ 2;
t82 = t89 ^ 2;
t77 = t82 * t95;
t71 = -t78 - t94;
t70 = -t77 - t94;
t66 = t77 + t78;
t65 = -t90 * qJDD(3) - t93 * t95;
t64 = t93 * qJDD(3) - t90 * t95;
t63 = (t82 + t83) * qJDD(3);
t59 = 0.2e1 * t103 + t107;
t48 = -t52 + t79;
t47 = t51 - t79;
t44 = -t52 - t79;
t43 = -t89 * t70 - t112;
t42 = t92 * t71 - t116;
t40 = t52 - t51;
t36 = -t79 - t51;
t35 = -t51 - t52;
t33 = -t55 * qJD(5) - t102;
t28 = -t88 * t44 - t113;
t27 = t91 * t44 - t117;
t26 = t50 + t34;
t21 = (qJD(5) + t81) * t55 + t102;
t17 = t91 * t36 - t126;
t16 = t88 * t36 + t125;
t12 = -t89 * t27 + t92 * t28;
t11 = t88 * t26 + t91 * t99;
t10 = -t91 * t26 + t88 * t99;
t8 = -t89 * t16 + t92 * t17;
t4 = -t89 * t10 + t92 * t11;
t1 = t92 * t3 - t121;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t46 + t86 * t97, 0, 0, 0, 0, 0, 0, t86 * t64 + t85 * t65, -t85 * t64 + t86 * t65, 0, t85 * (-t90 * t31 + t93 * t32) + t86 * (t93 * t31 + t90 * t32), 0, 0, 0, 0, 0, 0, t85 * (t93 * t42 - t90 * t61) + t86 * (t90 * t42 + t93 * t61), t85 * (t93 * t43 + t90 * t59) + t86 * (t90 * t43 - t93 * t59), t85 * (t93 * t63 - t90 * t66) + t86 * (t90 * t63 + t93 * t66), t85 * (t90 * t29 + t93 * t9) + t86 * (-t93 * t29 + t90 * t9), 0, 0, 0, 0, 0, 0, t85 * (t90 * t21 + t93 * t8) + t86 * (-t93 * t21 + t90 * t8), t85 * (t93 * t12 + t123 * t90) + t86 * (t90 * t12 - t123 * t93), t85 * (t90 * t35 + t93 * t4) + t86 * (-t93 * t35 + t90 * t4), t85 * (t93 * t1 + t90 * t15) + t86 * (t90 * t1 - t93 * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, t106 * t92 + t89 * t71, -t89 * t68 + t92 * t70, 0, -t92 * t19 + t89 * t20, 0, 0, 0, 0, 0, 0, t92 * t16 + t89 * t17, t92 * t27 + t89 * t28, t92 * t10 + t89 * t11, t92 * t2 + t89 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t31, -t32, 0, 0, (t60 + t103) * t89, t92 * t59 + t89 * t61, t116 + t92 * (-t77 + t94), t61 * t92, t89 * (t78 - t94) + t112, 0, pkin(3) * t61 + pkin(6) * t42 - t92 * t29, -pkin(3) * t59 + pkin(6) * t43 + t89 * t29, pkin(3) * t66 + pkin(6) * t63 + t9, -pkin(3) * t29 + pkin(6) * t9, t89 * (-t55 * t120 + t91 * t34) + t92 * (t55 * t119 + t88 * t34), t89 * (-t123 * t88 - t91 * t21) + t92 * (t123 * t91 - t88 * t21), t89 * (-t88 * t48 + t125) + t92 * (t91 * t48 + t126), t89 * (t53 * t119 - t88 * t33) + t92 * (t53 * t120 + t91 * t33), t89 * (t91 * t47 - t117) + t92 * (t88 * t47 + t113), (t89 * (-t53 * t91 + t55 * t88) + t92 * (-t53 * t88 - t55 * t91)) * t81, t89 * (-pkin(7) * t16 + t118) + t92 * (-pkin(4) * t21 + pkin(7) * t17 - t114) - pkin(3) * t21 + pkin(6) * t8, t89 * (-pkin(7) * t27 + t114) + t92 * (-pkin(4) * t123 + pkin(7) * t28 + t118) - pkin(3) * t123 + pkin(6) * t12, t89 * (-pkin(7) * t10 - t2) + t92 * (-pkin(4) * t35 + pkin(7) * t11 + t3) - pkin(3) * t35 + pkin(6) * t4, -pkin(7) * t121 + t92 * (-pkin(4) * t15 + pkin(7) * t3) - pkin(3) * t15 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t77 - t78, t107, t72, t76, qJDD(4), -t19, -t20, 0, 0, t41, t40, t26, -t41, t99, t80, pkin(4) * t16 - t6, -t115 - t88 * t122 + (-t106 * t88 + t27) * pkin(4), pkin(4) * t10, pkin(4) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, t26, -t41, t99, t80, -t6, -t7, 0, 0;];
tauJ_reg = t5;
