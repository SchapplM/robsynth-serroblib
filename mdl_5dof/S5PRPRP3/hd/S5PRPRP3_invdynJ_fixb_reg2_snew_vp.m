% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRP3
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:44
% DurationCPUTime: 0.52s
% Computational Cost: add. (1057->125), mult. (1990->170), div. (0->0), fcn. (1188->8), ass. (0->87)
t102 = (qJD(2) * qJD(5));
t115 = 2 * t102;
t86 = sin(qJ(4));
t88 = cos(qJ(4));
t91 = qJD(2) ^ 2;
t69 = t86 * t91 * t88;
t63 = qJDD(4) + t69;
t114 = pkin(4) * t63;
t82 = sin(pkin(7));
t84 = cos(pkin(7));
t61 = -t84 * g(1) - t82 * g(2);
t78 = -g(3) + qJDD(1);
t87 = sin(qJ(2));
t89 = cos(qJ(2));
t42 = t89 * t61 + t87 * t78;
t40 = -t91 * pkin(2) + t42;
t81 = sin(pkin(8));
t83 = cos(pkin(8));
t41 = -t87 * t61 + t89 * t78;
t92 = qJDD(2) * pkin(2) + t41;
t19 = t83 * t40 + t81 * t92;
t17 = -t91 * pkin(3) + qJDD(2) * pkin(6) + t19;
t55 = -t82 * g(1) + t84 * g(2) + qJDD(3);
t13 = t88 * t17 + t86 * t55;
t75 = t88 * qJDD(2);
t103 = qJD(2) * qJD(4);
t97 = t86 * t103;
t52 = t75 - t97;
t106 = qJD(2) * t86;
t62 = qJD(4) * pkin(4) - qJ(5) * t106;
t93 = t52 * qJ(5) - qJD(4) * t62 + t88 * t115 + t13;
t76 = t86 ^ 2;
t113 = t76 * t91;
t77 = t88 ^ 2;
t112 = t77 * t91;
t111 = t86 * t17;
t110 = t86 * t63;
t64 = qJDD(4) - t69;
t109 = t88 * t64;
t108 = t76 + t77;
t107 = qJ(5) * t86;
t105 = t86 * qJDD(2);
t104 = qJ(5) * qJDD(2);
t90 = qJD(4) ^ 2;
t67 = -t90 - t112;
t38 = t88 * t67 - t110;
t53 = t75 - 0.2e1 * t97;
t22 = t81 * t38 + t83 * t53;
t101 = pkin(2) * t22 + pkin(3) * t53 + pkin(6) * t38;
t66 = -t90 - t113;
t39 = -t86 * t66 - t109;
t96 = t88 * t103;
t50 = 0.2e1 * t96 + t105;
t23 = t81 * t39 - t83 * t50;
t100 = pkin(2) * t23 - pkin(3) * t50 + pkin(6) * t39;
t58 = t108 * qJDD(2);
t59 = t108 * t91;
t26 = t81 * t58 + t83 * t59;
t99 = pkin(2) * t26 + pkin(3) * t59 + pkin(6) * t58;
t44 = t88 * t55;
t12 = -t44 + t111;
t4 = t86 * t12 + t88 * t13;
t18 = -t81 * t40 + t83 * t92;
t16 = -qJDD(2) * pkin(3) - t91 * pkin(6) - t18;
t51 = t96 + t105;
t94 = -t44 + (t51 - t96) * qJ(5) - t114;
t5 = -0.2e1 * t86 * t102 - t111 - t94;
t8 = -t52 * pkin(4) - qJ(5) * t112 + t62 * t106 + qJDD(5) + t16;
t60 = (t76 - t77) * t91;
t57 = -t81 * qJDD(2) - t83 * t91;
t56 = t83 * qJDD(2) - t81 * t91;
t37 = -t86 * t64 + t88 * t66;
t36 = t110 + t88 * (t90 - t113);
t35 = t88 * t63 + t86 * t67;
t34 = t86 * (-t90 + t112) + t109;
t33 = (t51 + t96) * t86;
t32 = (t52 - t97) * t88;
t24 = t88 * t50 + t86 * t53;
t14 = t87 * (t83 * t58 - t81 * t59) + t89 * t26;
t10 = t87 * (t83 * t39 + t81 * t50) + t89 * t23;
t9 = t87 * (t83 * t38 - t81 * t53) + t89 * t22;
t7 = -pkin(4) * t112 + t93;
t6 = t83 * t18 + t81 * t19;
t3 = -t83 * t16 + t81 * t4;
t2 = -t86 * t5 + t88 * t7;
t1 = t81 * t2 - t83 * t8;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, 0, 0, 0, 0, 0, t89 * qJDD(2) - t87 * t91, -t87 * qJDD(2) - t89 * t91, 0, t89 * t41 + t87 * t42, 0, 0, 0, 0, 0, 0, t89 * t56 + t87 * t57, -t87 * t56 + t89 * t57, 0, t87 * (-t81 * t18 + t83 * t19) + t89 * t6, 0, 0, 0, 0, 0, 0, t9, t10, t14, t87 * (t81 * t16 + t83 * t4) + t89 * t3, 0, 0, 0, 0, 0, 0, t9, t10, t14, t87 * (t83 * t2 + t81 * t8) + t89 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t41, -t42, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t56 + t18, pkin(2) * t57 - t19, 0, pkin(2) * t6, t33, t24, t36, t32, t34, 0, -t88 * t16 + t101, t86 * t16 + t100, t4 + t99, pkin(2) * t3 - pkin(3) * t16 + pkin(6) * t4, t33, t24, t36, t32, t34, 0, -t63 * t107 + t88 * (pkin(4) * t53 + qJ(5) * t67 - t8) + t101, t86 * (-qJ(5) * t66 + t8) + t88 * (-pkin(4) * t50 - qJ(5) * t64) + t100, t88 * (t88 * t104 + (t59 - t112) * pkin(4) + t93) + ((t115 + t17 + t104) * t86 + t94) * t86 + t99, -t5 * t107 + t88 * (-pkin(4) * t8 + qJ(5) * t7) - pkin(3) * t8 + pkin(6) * t2 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, t35, t37, 0, -t88 * t12 + t86 * t13, 0, 0, 0, 0, 0, 0, t35, t37, 0, t88 * t5 + t86 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t60, t105, t69, t75, qJDD(4), -t12, -t13, 0, 0, -t69, t60, t105, t69, t75, qJDD(4), t5 + t114, (t66 + t112) * pkin(4) - t93, -pkin(4) * t105, pkin(4) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t50, -t59, t8;];
tauJ_reg = t11;
