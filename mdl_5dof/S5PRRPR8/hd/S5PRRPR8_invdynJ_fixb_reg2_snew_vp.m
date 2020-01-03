% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:44
% DurationCPUTime: 0.60s
% Computational Cost: add. (2270->110), mult. (3076->171), div. (0->0), fcn. (2072->10), ass. (0->88)
t101 = cos(qJ(3));
t89 = qJD(2) + qJD(3);
t87 = t89 ^ 2;
t88 = qJDD(2) + qJDD(3);
t93 = sin(pkin(9));
t95 = cos(pkin(9));
t107 = -t95 * t87 - t93 * t88;
t67 = t93 * t87 - t95 * t88;
t98 = sin(qJ(3));
t40 = t101 * t107 + t98 * t67;
t120 = t101 * t67 - t107 * t98;
t100 = cos(qJ(5));
t94 = sin(pkin(8));
t96 = cos(pkin(8));
t105 = -t94 * g(1) + t96 * g(2) + qJDD(4);
t102 = cos(qJ(2));
t76 = -t96 * g(1) - t94 * g(2);
t92 = -g(3) + qJDD(1);
t99 = sin(qJ(2));
t55 = t102 * t92 - t99 * t76;
t53 = qJDD(2) * pkin(2) + t55;
t104 = qJD(2) ^ 2;
t56 = t102 * t76 + t99 * t92;
t54 = -t104 * pkin(2) + t56;
t32 = t101 * t53 - t98 * t54;
t106 = t88 * pkin(3) + t32;
t33 = t101 * t54 + t98 * t53;
t29 = -t87 * pkin(3) + t33;
t21 = t93 * t106 + t95 * t29;
t19 = -t87 * pkin(4) + t88 * pkin(7) + t21;
t97 = sin(qJ(5));
t13 = -t100 * t105 + t97 * t19;
t14 = t100 * t19 + t97 * t105;
t7 = t100 * t14 + t97 * t13;
t79 = t97 * t87 * t100;
t74 = qJDD(5) + t79;
t117 = t97 * t74;
t116 = t97 * t88;
t75 = qJDD(5) - t79;
t115 = t100 * t75;
t114 = qJD(5) * t89;
t20 = t95 * t106 - t93 * t29;
t18 = -t88 * pkin(4) - t87 * pkin(7) - t20;
t4 = -t95 * t18 + t93 * t7;
t113 = pkin(3) * t4 - pkin(4) * t18 + pkin(7) * t7;
t112 = pkin(3) * t107 - t21;
t103 = qJD(5) ^ 2;
t90 = t97 ^ 2;
t82 = t90 * t87;
t77 = -t82 - t103;
t51 = -t97 * t77 - t115;
t62 = 0.2e1 * t100 * t114 + t116;
t35 = t93 * t51 - t95 * t62;
t111 = pkin(3) * t35 - pkin(4) * t62 + pkin(7) * t51 + t97 * t18;
t91 = t100 ^ 2;
t83 = t91 * t87;
t78 = -t83 - t103;
t50 = t100 * t78 - t117;
t81 = t100 * t88;
t63 = -0.2e1 * t97 * t114 + t81;
t34 = t93 * t50 + t95 * t63;
t110 = pkin(3) * t34 + pkin(4) * t63 + pkin(7) * t50 - t100 * t18;
t109 = -pkin(3) * t67 + t20;
t69 = (t90 + t91) * t88;
t73 = t82 + t83;
t42 = t93 * t69 + t95 * t73;
t108 = pkin(3) * t42 + pkin(4) * t73 + pkin(7) * t69 + t7;
t71 = t101 * t88 - t98 * t87;
t70 = -t101 * t87 - t98 * t88;
t49 = t117 + t100 * (-t82 + t103);
t48 = t97 * (t83 - t103) + t115;
t45 = t62 * t97;
t44 = t63 * t100;
t43 = t95 * t69 - t93 * t73;
t38 = t100 * t62 + t97 * t63;
t37 = t95 * t51 + t93 * t62;
t36 = t95 * t50 - t93 * t63;
t25 = t101 * t42 + t98 * t43;
t24 = t101 * t35 + t98 * t37;
t23 = t101 * t34 + t98 * t36;
t22 = t101 * t32 + t98 * t33;
t10 = -t93 * t20 + t95 * t21;
t9 = t95 * t20 + t93 * t21;
t8 = pkin(3) * t9;
t5 = t93 * t18 + t95 * t7;
t2 = t98 * t10 + t101 * t9;
t1 = t101 * t4 + t98 * t5;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, t102 * qJDD(2) - t99 * t104, -t99 * qJDD(2) - t102 * t104, 0, t102 * t55 + t99 * t56, 0, 0, 0, 0, 0, 0, t102 * t71 + t99 * t70, t102 * t70 - t99 * t71, 0, t99 * (t101 * t33 - t98 * t32) + t102 * t22, 0, 0, 0, 0, 0, 0, -t102 * t120 + t99 * t40, t102 * t40 + t99 * t120, 0, t99 * (t101 * t10 - t98 * t9) + t102 * t2, 0, 0, 0, 0, 0, 0, t99 * (t101 * t36 - t98 * t34) + t102 * t23, t99 * (t101 * t37 - t98 * t35) + t102 * t24, t99 * (t101 * t43 - t98 * t42) + t102 * t25, t99 * (t101 * t5 - t98 * t4) + t102 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t55, -t56, 0, 0, 0, 0, 0, 0, 0, t88, pkin(2) * t71 + t32, pkin(2) * t70 - t33, 0, pkin(2) * t22, 0, 0, 0, 0, 0, t88, -pkin(2) * t120 + t109, pkin(2) * t40 + t112, 0, pkin(2) * t2 + t8, t45, t38, t49, t44, t48, 0, pkin(2) * t23 + t110, pkin(2) * t24 + t111, pkin(2) * t25 + t108, pkin(2) * t1 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t32, -t33, 0, 0, 0, 0, 0, 0, 0, t88, t109, t112, 0, t8, t45, t38, t49, t44, t48, 0, t110, t111, t108, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, 0, 0, 0, t100 * t74 + t97 * t78, t100 * t77 - t97 * t75, 0, -t100 * t13 + t97 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t82 - t83, t116, t79, t81, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t3;
