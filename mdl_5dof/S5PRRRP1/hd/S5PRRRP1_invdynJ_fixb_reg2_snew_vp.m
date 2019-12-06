% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:12
% EndTime: 2019-12-05 16:40:15
% DurationCPUTime: 0.51s
% Computational Cost: add. (1520->112), mult. (2260->137), div. (0->0), fcn. (1365->8), ass. (0->88)
t81 = qJD(2) + qJD(3);
t110 = (qJD(5) * t81);
t125 = 2 * t110;
t79 = t81 ^ 2;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t67 = t88 * t79 * t91;
t59 = qJDD(4) + t67;
t124 = pkin(4) * t59;
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t62 = t86 * g(1) - t87 * g(2);
t63 = -t87 * g(1) - t86 * g(2);
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t103 = t93 * t62 - t90 * t63;
t36 = qJDD(2) * pkin(2) + t103;
t100 = -t90 * t62 - t93 * t63;
t37 = -qJD(2) ^ 2 * pkin(2) - t100;
t89 = sin(qJ(3));
t92 = cos(qJ(3));
t25 = t89 * t36 + t92 * t37;
t80 = qJDD(2) + qJDD(3);
t23 = -t79 * pkin(3) + t80 * pkin(7) + t25;
t120 = t88 * t23;
t84 = -g(3) + qJDD(1);
t76 = t91 * t84;
t16 = -t76 + t120;
t17 = t91 * t23 + t88 * t84;
t5 = t88 * t16 + t91 * t17;
t111 = qJD(4) * t81;
t106 = t88 * t111;
t74 = t91 * t80;
t51 = t74 - t106;
t112 = qJ(5) * t88;
t58 = qJD(4) * pkin(4) - t81 * t112;
t95 = t51 * qJ(5) - qJD(4) * t58 + t91 * t125 + t17;
t24 = t92 * t36 - t89 * t37;
t22 = -t80 * pkin(3) - t79 * pkin(7) - t24;
t123 = -pkin(3) * t22 + pkin(7) * t5;
t82 = t88 ^ 2;
t122 = t82 * t79;
t83 = t91 ^ 2;
t121 = t83 * t79;
t119 = t88 * t59;
t73 = t88 * t80;
t60 = qJDD(4) - t67;
t118 = t91 * t60;
t94 = qJD(4) ^ 2;
t65 = -t94 - t121;
t42 = t91 * t65 - t119;
t52 = t74 - 0.2e1 * t106;
t117 = pkin(3) * t52 + pkin(7) * t42;
t64 = -t94 - t122;
t43 = -t88 * t64 - t118;
t105 = t91 * t111;
t49 = t73 + 0.2e1 * t105;
t116 = -pkin(3) * t49 + pkin(7) * t43;
t114 = t82 + t83;
t54 = t114 * t80;
t55 = t114 * t79;
t115 = pkin(3) * t55 + pkin(7) * t54;
t113 = qJ(5) * t80;
t109 = t88 * t22 + t116;
t108 = -t91 * t22 + t117;
t50 = t73 + t105;
t96 = -t76 + (-t105 + t50) * qJ(5) - t124;
t104 = t88 * ((t125 + t23 + t113) * t88 + t96) + t91 * (t91 * t113 + (t55 - t121) * pkin(4) + t95) + t115;
t13 = t88 * t81 * t58 - t51 * pkin(4) - qJ(5) * t121 + qJDD(5) + t22;
t102 = t88 * (-qJ(5) * t64 + t13) + t91 * (-pkin(4) * t49 - qJ(5) * t60) + t116;
t101 = t115 + t5;
t10 = -0.2e1 * t88 * t110 - t120 - t96;
t11 = -pkin(4) * t121 + t95;
t2 = -t88 * t10 + t91 * t11;
t98 = pkin(7) * t2 - t10 * t112 - pkin(3) * t13 + t91 * (-pkin(4) * t13 + qJ(5) * t11);
t97 = -t59 * t112 + t117 + t91 * (pkin(4) * t52 + qJ(5) * t65 - t13);
t56 = (t82 - t83) * t79;
t41 = -t88 * t60 + t91 * t64;
t40 = t119 + t91 * (t94 - t122);
t39 = t91 * t59 + t88 * t65;
t38 = t88 * (-t94 + t121) + t118;
t33 = (t50 + t105) * t88;
t32 = (t51 - t106) * t91;
t29 = pkin(2) * (t89 * t54 + t92 * t55);
t28 = t91 * t49 + t88 * t52;
t27 = pkin(2) * (t89 * t43 - t92 * t49);
t26 = pkin(2) * (t89 * t42 + t92 * t52);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, t39, t41, 0, -t91 * t16 + t88 * t17, 0, 0, 0, 0, 0, 0, t39, t41, 0, t91 * t10 + t88 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t103, t100, 0, 0, 0, 0, 0, 0, 0, t80, pkin(2) * (-t89 * t79 + t92 * t80) + t24, pkin(2) * (-t92 * t79 - t89 * t80) - t25, 0, pkin(2) * (t92 * t24 + t89 * t25), t33, t28, t40, t32, t38, 0, t26 + t108, t27 + t109, t29 + t101, pkin(2) * (-t92 * t22 + t89 * t5) + t123, t33, t28, t40, t32, t38, 0, t26 + t97, t27 + t102, t29 + t104, pkin(2) * (-t92 * t13 + t89 * t2) + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t24, -t25, 0, 0, t33, t28, t40, t32, t38, 0, t108, t109, t101, t123, t33, t28, t40, t32, t38, 0, t97, t102, t104, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t56, t73, t67, t74, qJDD(4), -t16, -t17, 0, 0, -t67, t56, t73, t67, t74, qJDD(4), t10 + t124, (t64 + t121) * pkin(4) - t95, -pkin(4) * t73, pkin(4) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t49, -t55, t13;];
tauJ_reg = t1;
