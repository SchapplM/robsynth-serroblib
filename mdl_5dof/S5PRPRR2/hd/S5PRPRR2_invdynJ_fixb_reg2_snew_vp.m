% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:16
% DurationCPUTime: 0.56s
% Computational Cost: add. (1953->110), mult. (2819->171), div. (0->0), fcn. (1902->10), ass. (0->85)
t83 = qJD(2) + qJD(4);
t81 = t83 ^ 2;
t82 = qJDD(2) + qJDD(4);
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t100 = -t95 * t81 - t92 * t82;
t61 = t92 * t81 - t95 * t82;
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t33 = t100 * t89 + t87 * t61;
t111 = -t100 * t87 + t89 * t61;
t88 = sin(pkin(8));
t90 = cos(pkin(8));
t70 = -t90 * g(1) - t88 * g(2);
t86 = -g(3) + qJDD(1);
t93 = sin(qJ(2));
t96 = cos(qJ(2));
t50 = -t93 * t70 + t96 * t86;
t48 = qJDD(2) * pkin(2) + t50;
t51 = t96 * t70 + t93 * t86;
t98 = qJD(2) ^ 2;
t49 = -t98 * pkin(2) + t51;
t27 = t87 * t48 + t89 * t49;
t25 = -t98 * pkin(3) + t27;
t26 = t89 * t48 - t87 * t49;
t99 = qJDD(2) * pkin(3) + t26;
t19 = t95 * t25 + t92 * t99;
t17 = -t81 * pkin(4) + t82 * pkin(7) + t19;
t67 = -t88 * g(1) + t90 * g(2) + qJDD(3);
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t11 = t91 * t17 - t94 * t67;
t12 = t94 * t17 + t91 * t67;
t6 = t91 * t11 + t94 * t12;
t18 = -t92 * t25 + t95 * t99;
t16 = -t82 * pkin(4) - t81 * pkin(7) - t18;
t108 = -pkin(4) * t16 + pkin(7) * t6;
t73 = t91 * t81 * t94;
t65 = qJDD(5) + t73;
t107 = t91 * t65;
t106 = t91 * t82;
t66 = qJDD(5) - t73;
t105 = t94 * t66;
t104 = qJD(5) * t83;
t84 = t91 ^ 2;
t76 = t84 * t81;
t97 = qJD(5) ^ 2;
t71 = -t76 - t97;
t44 = -t91 * t71 - t105;
t55 = 0.2e1 * t94 * t104 + t106;
t103 = -pkin(4) * t55 + pkin(7) * t44 + t91 * t16;
t85 = t94 ^ 2;
t77 = t85 * t81;
t72 = -t77 - t97;
t43 = t94 * t72 - t107;
t75 = t94 * t82;
t56 = -0.2e1 * t91 * t104 + t75;
t102 = pkin(4) * t56 + pkin(7) * t43 - t94 * t16;
t58 = (t84 + t85) * t82;
t64 = t76 + t77;
t101 = pkin(4) * t64 + pkin(7) * t58 + t6;
t69 = -t87 * qJDD(2) - t89 * t98;
t68 = t89 * qJDD(2) - t87 * t98;
t42 = t107 + t94 * (-t76 + t97);
t41 = t91 * (t77 - t97) + t105;
t38 = t55 * t91;
t37 = t56 * t94;
t36 = t95 * t58 - t92 * t64;
t35 = t92 * t58 + t95 * t64;
t32 = t94 * t55 + t91 * t56;
t31 = t95 * t44 + t92 * t55;
t30 = t95 * t43 - t92 * t56;
t29 = t92 * t44 - t95 * t55;
t28 = t92 * t43 + t95 * t56;
t23 = t89 * t35 + t87 * t36;
t22 = t89 * t29 + t87 * t31;
t21 = t89 * t28 + t87 * t30;
t20 = t89 * t26 + t87 * t27;
t8 = -t92 * t18 + t95 * t19;
t7 = t95 * t18 + t92 * t19;
t4 = t92 * t16 + t95 * t6;
t3 = -t95 * t16 + t92 * t6;
t2 = t89 * t7 + t87 * t8;
t1 = t89 * t3 + t87 * t4;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, t96 * qJDD(2) - t93 * t98, -t93 * qJDD(2) - t96 * t98, 0, t96 * t50 + t93 * t51, 0, 0, 0, 0, 0, 0, t96 * t68 + t93 * t69, -t93 * t68 + t96 * t69, 0, t93 * (-t87 * t26 + t89 * t27) + t96 * t20, 0, 0, 0, 0, 0, 0, -t111 * t96 + t93 * t33, t93 * t111 + t96 * t33, 0, t93 * (-t87 * t7 + t89 * t8) + t96 * t2, 0, 0, 0, 0, 0, 0, t93 * (-t87 * t28 + t89 * t30) + t96 * t21, t93 * (-t87 * t29 + t89 * t31) + t96 * t22, t93 * (-t87 * t35 + t89 * t36) + t96 * t23, t93 * (-t87 * t3 + t89 * t4) + t96 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t50, -t51, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t68 + t26, pkin(2) * t69 - t27, 0, pkin(2) * t20, 0, 0, 0, 0, 0, t82, -pkin(2) * t111 - pkin(3) * t61 + t18, pkin(2) * t33 + pkin(3) * t100 - t19, 0, pkin(2) * t2 + pkin(3) * t7, t38, t32, t42, t37, t41, 0, pkin(2) * t21 + pkin(3) * t28 + t102, pkin(2) * t22 + pkin(3) * t29 + t103, pkin(2) * t23 + pkin(3) * t35 + t101, pkin(2) * t1 + pkin(3) * t3 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, t94 * t65 + t91 * t72, -t91 * t66 + t94 * t71, 0, -t94 * t11 + t91 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t18, -t19, 0, 0, t38, t32, t42, t37, t41, 0, t102, t103, t101, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t76 - t77, t106, t73, t75, qJDD(5), -t11, -t12, 0, 0;];
tauJ_reg = t5;
