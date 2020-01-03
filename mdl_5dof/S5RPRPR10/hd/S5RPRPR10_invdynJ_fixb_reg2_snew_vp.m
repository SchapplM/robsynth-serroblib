% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:07
% DurationCPUTime: 0.62s
% Computational Cost: add. (2932->128), mult. (4144->159), div. (0->0), fcn. (1698->8), ass. (0->92)
t73 = (-qJD(1) + qJD(3));
t71 = t73 ^ 2;
t72 = qJDD(1) - qJDD(3);
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t49 = t77 * t71 + t78 * t72;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t96 = -t78 * t71 + t77 * t72;
t30 = t80 * t49 + t83 * t96;
t114 = t83 * t49 - t80 * t96;
t110 = pkin(1) + pkin(2);
t86 = qJD(1) ^ 2;
t74 = qJDD(1) * qJ(2);
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t99 = t84 * g(1) + t81 * g(2);
t95 = (2 * qJD(2) * qJD(1)) - t99;
t92 = t74 + t95;
t42 = -t110 * t86 + t92;
t106 = t81 * g(1) - t84 * g(2);
t102 = -qJDD(2) + t106;
t93 = -t86 * qJ(2) - t102;
t43 = -t110 * qJDD(1) + t93;
t28 = t83 * t42 + t80 * t43;
t22 = -t71 * pkin(3) + t28;
t97 = t80 * t42 - t83 * t43;
t90 = -t72 * pkin(3) - t97;
t101 = t77 * t22 - t78 * t90;
t14 = t78 * t22 + t77 * t90;
t6 = -t101 * t78 + t77 * t14;
t111 = pkin(3) * t6;
t79 = sin(qJ(5));
t82 = cos(qJ(5));
t59 = t82 * t71 * t79;
t55 = qJDD(5) + t59;
t109 = t79 * t55;
t108 = t79 * t72;
t56 = qJDD(5) - t59;
t107 = t82 * t56;
t105 = qJDD(1) * pkin(1);
t104 = g(3) + qJDD(4);
t103 = 2 * qJD(5) * t73;
t12 = -t71 * pkin(4) - t72 * pkin(7) + t14;
t10 = t79 * t104 + t82 * t12;
t9 = -t82 * t104 + t79 * t12;
t5 = t82 * t10 + t79 * t9;
t100 = pkin(3) * t96 - t14;
t61 = t82 * t72;
t98 = t79 * t103 + t61;
t52 = -t83 * t71 + t80 * t72;
t53 = t80 * t71 + t83 * t72;
t94 = -pkin(3) * t49 - t101;
t45 = t82 * t103 - t108;
t11 = t72 * pkin(4) - t71 * pkin(7) + t101;
t3 = -t78 * t11 + t77 * t5;
t91 = pkin(3) * t3 - pkin(4) * t11 + pkin(7) * t5;
t76 = t82 ^ 2;
t64 = t76 * t71;
t85 = qJD(5) ^ 2;
t58 = -t64 - t85;
t38 = t82 * t58 - t109;
t23 = t77 * t38 - t78 * t98;
t89 = pkin(3) * t23 - pkin(4) * t98 + pkin(7) * t38 - t82 * t11;
t75 = t79 ^ 2;
t63 = t75 * t71;
t57 = -t63 - t85;
t39 = -t79 * t57 - t107;
t24 = t77 * t39 - t78 * t45;
t88 = pkin(3) * t24 - pkin(4) * t45 + pkin(7) * t39 + t79 * t11;
t51 = (-t75 - t76) * t72;
t54 = t63 + t64;
t32 = t77 * t51 + t78 * t54;
t87 = pkin(3) * t32 + pkin(4) * t54 + pkin(7) * t51 + t5;
t44 = -t93 + t105;
t37 = t109 + t82 * (-t63 + t85);
t36 = t79 * (t64 - t85) + t107;
t35 = t45 * t79;
t34 = t98 * t82;
t33 = t78 * t51 - t77 * t54;
t29 = t82 * t45 - t79 * t98;
t26 = t78 * t39 + t77 * t45;
t25 = t78 * t38 + t77 * t98;
t18 = t83 * t32 + t80 * t33;
t17 = t80 * t28 - t83 * t97;
t16 = t83 * t24 + t80 * t26;
t15 = t83 * t23 + t80 * t25;
t7 = t101 * t77 + t78 * t14;
t4 = t77 * t11 + t78 * t5;
t2 = t83 * t6 + t80 * t7;
t1 = t83 * t3 + t80 * t4;
t8 = [0, 0, 0, 0, 0, qJDD(1), t106, t99, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t102 + 0.2e1 * t105, 0, 0.2e1 * t74 + t95, qJ(2) * (-t86 * pkin(1) + t92) + pkin(1) * t44, 0, 0, 0, 0, 0, t72, qJ(2) * t52 + t110 * t53 + t97, qJ(2) * t53 - t110 * t52 + t28, 0, qJ(2) * (t83 * t28 + t80 * t97) - t110 * t17, 0, 0, 0, 0, 0, t72, qJ(2) * t30 + t110 * t114 - t94, qJ(2) * t114 - t110 * t30 - t100, 0, qJ(2) * (-t80 * t6 + t83 * t7) - t111 - t110 * t2, -t35, -t29, -t37, t34, -t36, 0, qJ(2) * (-t80 * t23 + t83 * t25) - t110 * t15 - t89, qJ(2) * (-t80 * t24 + t83 * t26) - t110 * t16 - t88, qJ(2) * (-t80 * t32 + t83 * t33) - t110 * t18 - t87, qJ(2) * (-t80 * t3 + t83 * t4) - t110 * t1 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t86, -t44, 0, 0, 0, 0, 0, 0, -t53, t52, 0, t17, 0, 0, 0, 0, 0, 0, -t114, t30, 0, t2, 0, 0, 0, 0, 0, 0, t15, t16, t18, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t97, -t28, 0, 0, 0, 0, 0, 0, 0, -t72, t94, t100, 0, t111, t35, t29, t37, -t34, t36, 0, t89, t88, t87, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, t82 * t55 + t79 * t58, -t79 * t56 + t82 * t57, 0, t79 * t10 - t82 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t63 - t64, -t108, t59, -t61, qJDD(5), -t9, -t10, 0, 0;];
tauJ_reg = t8;
