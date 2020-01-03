% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:32
% DurationCPUTime: 0.76s
% Computational Cost: add. (2103->154), mult. (2964->157), div. (0->0), fcn. (1211->6), ass. (0->101)
t140 = pkin(1) + pkin(2);
t80 = -qJD(1) + qJD(3);
t78 = t80 ^ 2;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t66 = t90 * t78 * t87;
t61 = qJDD(4) - t66;
t129 = t90 * t61;
t84 = t87 ^ 2;
t133 = t84 * t78;
t93 = qJD(4) ^ 2;
t62 = t93 + t133;
t38 = -t87 * t62 + t129;
t121 = qJD(4) * t80;
t118 = t90 * t121;
t79 = qJDD(1) - qJDD(3);
t130 = t87 * t79;
t52 = 0.2e1 * t118 - t130;
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t22 = t88 * t38 + t91 * t52;
t147 = t140 * t22;
t146 = pkin(7) * t38;
t145 = qJ(5) * t52;
t60 = qJDD(4) + t66;
t131 = t87 * t60;
t85 = t90 ^ 2;
t132 = t85 * t78;
t64 = -t93 - t132;
t37 = t90 * t64 - t131;
t119 = t87 * t121;
t128 = t90 * t79;
t53 = -0.2e1 * t119 - t128;
t21 = t88 * t37 + t91 * t53;
t144 = qJ(2) * (t91 * t37 - t88 * t53) - t140 * t21;
t124 = t84 + t85;
t55 = t124 * t79;
t58 = t124 * t78;
t30 = -t88 * t55 + t91 * t58;
t143 = qJ(2) * (-t91 * t55 - t88 * t58) - t140 * t30;
t122 = t87 * qJ(5);
t136 = t90 * pkin(4);
t106 = -t122 - t136;
t135 = t106 * t78;
t89 = sin(qJ(1));
t92 = cos(qJ(1));
t108 = t92 * g(1) + t89 * g(2);
t104 = (2 * qJD(2) * qJD(1)) - t108;
t82 = qJDD(1) * qJ(2);
t101 = t104 + t82;
t142 = qJD(1) ^ 2;
t44 = -t140 * t142 + t101;
t116 = t89 * g(1) - t92 * g(2);
t107 = qJDD(2) - t116;
t98 = -t142 * qJ(2) + t107;
t94 = -t140 * qJDD(1) + t98;
t25 = t91 * t44 + t88 * t94;
t18 = -t78 * pkin(3) - t79 * pkin(7) + t25;
t77 = t90 * g(3);
t12 = -qJDD(4) * pkin(4) - t93 * qJ(5) + (t18 + t135) * t87 + qJDD(5) - t77;
t141 = 2 * qJD(5);
t134 = t80 * t87;
t102 = -t119 - t128;
t113 = t88 * t44 - t91 * t94;
t17 = t79 * pkin(3) - t78 * pkin(7) + t113;
t96 = -t102 * pkin(4) - t145 + t17;
t95 = t134 * t141 - t96;
t137 = t87 * (-pkin(4) * t119 + t145 + t95);
t15 = t87 * g(3) + t90 * t18;
t126 = pkin(3) * t53 + pkin(7) * t37;
t125 = pkin(3) * t58 - pkin(7) * t55;
t123 = qJ(2) * t88;
t120 = qJDD(1) * pkin(1);
t117 = pkin(3) + t136;
t115 = qJ(2) * t91 - pkin(7);
t14 = t87 * t18 - t77;
t4 = t87 * t14 + t90 * t15;
t111 = -t90 * t17 + t126;
t110 = -pkin(3) * t17 + pkin(7) * t4;
t109 = qJDD(4) * qJ(5) + (qJD(4) * t141) + t90 * t135 + t15;
t26 = t90 * t52 + t87 * t53;
t34 = t87 * (-t93 + t132) + t129;
t56 = -t91 * t78 + t88 * t79;
t57 = t88 * t78 + t91 * t79;
t105 = t4 + t125;
t103 = -pkin(3) + t106;
t100 = t90 * ((t58 - t93) * pkin(4) + t109) + t87 * (qJ(5) * t58 + t12) + t125;
t99 = pkin(3) * t52 - t87 * t17 + t146;
t97 = t53 * t122 + t90 * ((t53 - t119) * pkin(4) + t95) + t126;
t59 = (t84 - t85) * t78;
t49 = -t98 + t120;
t36 = t131 + t90 * (t93 - t133);
t32 = t52 * t87;
t31 = -t90 * t102 + t87 * t118;
t11 = -t93 * pkin(4) + t109;
t10 = -t113 * t91 + t88 * t25;
t7 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t134 + t96;
t3 = t90 * t11 + t87 * t12;
t2 = -t91 * t17 + t88 * t4;
t1 = t88 * t3 - t91 * t7;
t5 = [0, 0, 0, 0, 0, qJDD(1), t116, t108, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t107 + 0.2e1 * t120, 0, t104 + 0.2e1 * t82, pkin(1) * t49 + qJ(2) * (-t142 * pkin(1) + t101), 0, 0, 0, 0, 0, t79, qJ(2) * t56 + t140 * t57 + t113, qJ(2) * t57 - t140 * t56 + t25, 0, qJ(2) * (t113 * t88 + t91 * t25) - t140 * t10, -t32, -t26, -t36, t31, -t34, 0, -t111 + t144, qJ(2) * (-t38 * t91 + t88 * t52) + t147 + t99, -t105 + t143, qJ(2) * (t88 * t17 + t91 * t4) - t140 * t2 - t110, -t32, -t36, t26, 0, t34, t31, -t97 + t144, -t100 + t143, -t137 + t115 * t38 + (-t117 - t123) * t52 - t147, t115 * t3 - t140 * t1 + (-t103 + t123) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t142, -t49, 0, 0, 0, 0, 0, 0, -t57, t56, 0, t10, 0, 0, 0, 0, 0, 0, t21, -t22, t30, t2, 0, 0, 0, 0, 0, 0, t21, t30, t22, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t113, -t25, 0, 0, t32, t26, t36, -t31, t34, 0, t111, -t99, t105, t110, t32, t36, -t26, 0, -t34, -t31, t97, t100, t117 * t52 + t137 + t146, pkin(7) * t3 + t103 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t59, -t130, t66, -t128, qJDD(4), -t14, -t15, 0, 0, -t66, -t130, -t59, qJDD(4), t128, t66, pkin(4) * t60 + qJ(5) * t64 - t12, (pkin(4) * t87 - qJ(5) * t90) * t79, qJ(5) * t61 + (t62 - t93) * pkin(4) + t109, -pkin(4) * t12 + qJ(5) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t130, -t62, t12;];
tauJ_reg = t5;
