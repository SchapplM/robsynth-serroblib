% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR8
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:18
% DurationCPUTime: 0.68s
% Computational Cost: add. (1380->149), mult. (2781->193), div. (0->0), fcn. (1694->6), ass. (0->104)
t74 = sin(qJ(4));
t75 = sin(qJ(3));
t77 = cos(qJ(4));
t78 = cos(qJ(3));
t50 = (-t74 * t78 - t75 * t77) * qJD(1);
t100 = qJD(1) * t78;
t52 = -qJD(1) * t74 * t75 + t100 * t77;
t117 = t52 * t50;
t69 = qJDD(3) + qJDD(4);
t122 = t69 + t117;
t125 = t122 * t74;
t124 = t122 * t77;
t70 = qJD(3) + qJD(4);
t116 = t70 * t50;
t97 = qJD(1) * qJD(3);
t92 = t78 * t97;
t98 = t75 * qJDD(1);
t58 = -t92 - t98;
t65 = t78 * qJDD(1);
t93 = t75 * t97;
t59 = t65 - t93;
t27 = qJD(4) * t50 + t58 * t74 + t59 * t77;
t123 = t27 + t116;
t81 = qJD(1) ^ 2;
t103 = t78 * t81;
t120 = pkin(5) + pkin(1);
t101 = t81 * qJ(2);
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t94 = t76 * g(1) - t79 * g(2);
t89 = qJDD(2) - t94;
t86 = t89 - t101;
t44 = -qJDD(1) * t120 + t86;
t105 = t78 * t44;
t119 = t59 * pkin(6);
t121 = (-pkin(3) * t103 - pkin(6) * t97 + g(3)) * t75 + qJDD(3) * pkin(3) + t105 - t119;
t48 = t50 ^ 2;
t49 = t52 ^ 2;
t68 = t70 ^ 2;
t72 = t75 ^ 2;
t113 = t72 * t81;
t36 = t78 * g(3) - t44 * t75;
t87 = -qJD(3) * pkin(3) + pkin(6) * t100;
t24 = -pkin(3) * t113 + t58 * pkin(6) + qJD(3) * t87 - t36;
t10 = -t121 * t77 + t74 * t24;
t108 = t77 * t24;
t11 = t121 * t74 + t108;
t3 = -t10 * t77 + t11 * t74;
t118 = t78 * t3;
t115 = t70 * t74;
t114 = t70 * t77;
t73 = t78 ^ 2;
t112 = t73 * t81;
t71 = qJDD(1) * qJ(2);
t90 = t79 * g(1) + t76 * g(2);
t88 = -t71 + t90;
t96 = qJD(2) * qJD(1);
t83 = t88 - 0.2e1 * t96;
t25 = t58 * pkin(3) + t87 * t100 + (pkin(6) * t72 + t120) * t81 + t83;
t111 = t74 * t25;
t31 = -t117 + t69;
t110 = t74 * t31;
t95 = t75 * t103;
t109 = t75 * (qJDD(3) + t95);
t107 = t77 * t25;
t106 = t77 * t31;
t63 = qJDD(3) - t95;
t104 = t78 * t63;
t102 = t72 + t73;
t99 = qJDD(1) * pkin(1);
t4 = t10 * t74 + t11 * t77;
t91 = -t58 * t77 + t74 * t59;
t35 = t75 * g(3) + t105;
t22 = t78 * t35 - t75 * t36;
t84 = (-qJD(4) + t70) * t52 - t91;
t80 = qJD(3) ^ 2;
t66 = 0.2e1 * t96;
t61 = t102 * qJDD(1);
t60 = t65 - 0.2e1 * t93;
t57 = 0.2e1 * t92 + t98;
t47 = -t86 + t99;
t43 = -t49 + t68;
t42 = t48 - t68;
t41 = t120 * t81 + t83;
t39 = -t49 - t68;
t38 = -t109 + t78 * (-t80 - t112);
t37 = t75 * (-t80 - t113) + t104;
t33 = t49 - t48;
t29 = -t68 - t48;
t28 = -t48 - t49;
t26 = -qJD(4) * t52 - t91;
t21 = -t39 * t74 - t106;
t20 = t39 * t77 - t110;
t19 = -t116 + t27;
t14 = (qJD(4) + t70) * t52 + t91;
t13 = t29 * t77 - t125;
t12 = t29 * t74 + t124;
t8 = t20 * t78 + t21 * t75;
t7 = t19 * t74 + t77 * t84;
t6 = -t19 * t77 + t74 * t84;
t5 = t12 * t78 + t13 * t75;
t2 = t6 * t78 + t7 * t75;
t1 = t4 * t75 + t118;
t9 = [0, 0, 0, 0, 0, qJDD(1), t94, t90, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t89 - 0.2e1 * t99, t66 + 0.2e1 * t71 - t90, pkin(1) * t47 + qJ(2) * (-pkin(1) * t81 + t66 - t88), (t59 - t93) * t78, -t57 * t78 - t60 * t75, t104 - t75 * (t80 - t112), (-t58 + t92) * t75, t78 * (-t80 + t113) - t109, 0, qJ(2) * t57 - t120 * t37 - t75 * t41, qJ(2) * t60 - t120 * t38 - t78 * t41, -t101 * t102 + t120 * t61 - t22, -qJ(2) * t41 - t120 * t22, t78 * (-t115 * t52 + t27 * t77) - t75 * (t114 * t52 + t27 * t74), t78 * (-t123 * t74 - t14 * t77) - t75 * (t123 * t77 - t14 * t74), t78 * (-t43 * t74 + t124) - t75 * (t43 * t77 + t125), t78 * (-t114 * t50 - t26 * t74) - t75 * (-t115 * t50 + t26 * t77), t78 * (t42 * t77 - t110) - t75 * (t42 * t74 + t106), (t78 * (t50 * t77 + t52 * t74) - t75 * (t50 * t74 - t52 * t77)) * t70, t78 * (-pkin(6) * t12 - t111) - t75 * (-pkin(3) * t14 + pkin(6) * t13 + t107) + qJ(2) * t14 - t120 * t5, t78 * (-pkin(6) * t20 - t107) - t75 * (-pkin(3) * t123 + pkin(6) * t21 - t111) + qJ(2) * t123 - t120 * t8, t78 * (-pkin(6) * t6 - t3) - t75 * (-pkin(3) * t28 + pkin(6) * t7 + t4) + qJ(2) * t28 - t120 * t2, -pkin(6) * t118 - t75 * (pkin(3) * t25 + pkin(6) * t4) - qJ(2) * t25 - t120 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t81, -t47, 0, 0, 0, 0, 0, 0, t37, t38, -t61, t22, 0, 0, 0, 0, 0, 0, t5, t8, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, (-t72 + t73) * t81, t65, -t95, -t98, qJDD(3), t35, t36, 0, 0, -t117, t33, t19, t117, t84, t69, pkin(3) * t12 - t10, -t108 - t74 * (-pkin(6) * t93 - t119 + t35) + (-t63 * t74 + t20) * pkin(3), pkin(3) * t6, pkin(3) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t33, t19, t117, t84, t69, -t10, -t11, 0, 0;];
tauJ_reg = t9;
