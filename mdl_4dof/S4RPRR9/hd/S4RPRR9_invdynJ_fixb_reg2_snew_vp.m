% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR9
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:30
% DurationCPUTime: 0.60s
% Computational Cost: add. (1239->141), mult. (2425->180), div. (0->0), fcn. (1418->6), ass. (0->106)
t76 = cos(qJ(3));
t104 = qJD(1) * t76;
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t52 = -t75 * qJD(3) + t72 * t104;
t54 = t72 * qJD(3) + t75 * t104;
t38 = t54 * t52;
t73 = sin(qJ(3));
t100 = t73 * qJDD(1);
t99 = qJD(1) * qJD(3);
t93 = t76 * t99;
t56 = -t93 - t100;
t51 = qJDD(4) - t56;
t119 = -t38 + t51;
t121 = t119 * t72;
t120 = t119 * t75;
t117 = pkin(5) + pkin(1);
t79 = qJD(1) ^ 2;
t118 = t117 * t79;
t63 = t73 * qJD(1) + qJD(4);
t65 = t76 * qJDD(1);
t94 = t73 * t99;
t57 = t65 - t94;
t92 = -t75 * qJDD(3) + t72 * t57;
t19 = (qJD(4) - t63) * t54 + t92;
t49 = t52 ^ 2;
t50 = t54 ^ 2;
t62 = t63 ^ 2;
t116 = t63 * t72;
t115 = t63 * t75;
t70 = t73 ^ 2;
t114 = t70 * t79;
t71 = t76 ^ 2;
t113 = t71 * t79;
t105 = t79 * qJ(2);
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t95 = t74 * g(1) - t77 * g(2);
t88 = qJDD(2) - t95;
t80 = t88 - t105;
t45 = -t117 * qJDD(1) + t80;
t34 = t73 * g(3) + t76 * t45;
t78 = qJD(3) ^ 2;
t91 = pkin(3) * t73 - pkin(6) * t76;
t82 = t79 * t91;
t25 = qJDD(3) * pkin(3) + t78 * pkin(6) - t76 * t82 + t34;
t112 = t72 * t25;
t29 = t38 + t51;
t111 = t72 * t29;
t96 = t73 * t79 * t76;
t110 = t73 * (qJDD(3) + t96);
t109 = t75 * t25;
t108 = t75 * t29;
t107 = t76 * (qJDD(3) - t96);
t106 = t70 + t71;
t103 = qJDD(1) * pkin(1);
t101 = qJD(4) + t63;
t98 = qJD(2) * qJD(1);
t97 = t73 * t38;
t67 = 0.2e1 * t98;
t69 = qJDD(1) * qJ(2);
t90 = t77 * g(1) + t74 * g(2);
t84 = -t69 + t90;
t81 = t67 - t84;
t86 = -t57 + t94;
t87 = -t56 + t93;
t18 = t87 * pkin(3) + t86 * pkin(6) - t118 + t81;
t35 = t76 * g(3) - t73 * t45;
t26 = -t78 * pkin(3) + qJDD(3) * pkin(6) - t73 * t82 - t35;
t8 = -t75 * t18 + t72 * t26;
t9 = t72 * t18 + t75 * t26;
t3 = t72 * t8 + t75 * t9;
t89 = t72 * t9 - t75 * t8;
t17 = t76 * t34 - t73 * t35;
t85 = -t72 * qJDD(3) - t75 * t57;
t83 = qJ(2) + t91;
t32 = -t52 * qJD(4) - t85;
t59 = t106 * qJDD(1);
t58 = t65 - 0.2e1 * t94;
t55 = 0.2e1 * t93 + t100;
t47 = -t80 + t103;
t46 = t63 * t52;
t44 = -t50 + t62;
t43 = t49 - t62;
t42 = t84 - 0.2e1 * t98 + t118;
t40 = -t110 + t76 * (-t78 - t113);
t39 = t73 * (-t78 - t114) + t107;
t37 = t50 - t49;
t36 = -t50 - t62;
t33 = -t62 - t49;
t31 = -t54 * qJD(4) - t92;
t27 = t49 + t50;
t24 = t101 * t52 + t85;
t23 = t32 + t46;
t22 = t32 - t46;
t20 = -t101 * t54 - t92;
t15 = -t72 * t36 - t108;
t14 = t75 * t36 - t111;
t13 = t75 * t33 - t121;
t12 = t72 * t33 + t120;
t11 = -t19 * t75 + t72 * t23;
t6 = t73 * t15 + t76 * t24;
t5 = t73 * t13 + t76 * t20;
t4 = t73 * t11 + t76 * t27;
t1 = t76 * t25 + t73 * t3;
t2 = [0, 0, 0, 0, 0, qJDD(1), t95, t90, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t88 - 0.2e1 * t103, t67 + 0.2e1 * t69 - t90, pkin(1) * t47 + qJ(2) * (-t79 * pkin(1) + t81), -t86 * t76, -t76 * t55 - t73 * t58, t107 - t73 * (t78 - t113), t87 * t73, t76 * (-t78 + t114) - t110, 0, qJ(2) * t55 - t117 * t39 - t73 * t42, qJ(2) * t58 - t117 * t40 - t76 * t42, -t106 * t105 + t117 * t59 - t17, -qJ(2) * t42 - t117 * t17, t76 * (-t54 * t116 + t75 * t32) + t97, t76 * (t75 * t20 - t72 * t22) + t73 * t37, t76 * (-t72 * t44 + t120) + t73 * t23, t76 * (t52 * t115 - t72 * t31) - t97, t76 * (t75 * t43 - t111) - t73 * t19, t73 * t51 + t76 * (-t52 * t75 + t54 * t72) * t63, t76 * (-pkin(6) * t12 - t112) - t73 * (-pkin(3) * t12 + t8) + qJ(2) * t12 - t117 * t5, t76 * (-pkin(6) * t14 - t109) - t73 * (-pkin(3) * t14 + t9) + qJ(2) * t14 - t117 * t6, -t76 * t89 - t117 * t4 + t83 * (-t19 * t72 - t75 * t23), -t117 * t1 + t83 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t79, -t47, 0, 0, 0, 0, 0, 0, t39, t40, -t59, t17, 0, 0, 0, 0, 0, 0, t5, t6, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, (-t70 + t71) * t79, t65, -t96, -t100, qJDD(3), t34, t35, 0, 0, t54 * t115 + t72 * t32, t72 * t20 + t75 * t22, t75 * t44 + t121, t52 * t116 + t75 * t31, t72 * t43 + t108, (-t52 * t72 - t54 * t75) * t63, pkin(3) * t20 + pkin(6) * t13 + t109, pkin(3) * t24 + pkin(6) * t15 - t112, pkin(3) * t27 + pkin(6) * t11 + t3, pkin(3) * t25 + pkin(6) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t37, t23, -t38, -t19, t51, -t8, -t9, 0, 0;];
tauJ_reg = t2;
