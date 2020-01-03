% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:22
% EndTime: 2020-01-03 11:52:26
% DurationCPUTime: 0.89s
% Computational Cost: add. (3491->125), mult. (4975->187), div. (0->0), fcn. (2890->10), ass. (0->97)
t100 = sin(qJ(3));
t104 = cos(qJ(3));
t103 = cos(qJ(4));
t92 = qJD(1) + qJD(3);
t89 = qJD(4) + t92;
t87 = t89 ^ 2;
t91 = qJDD(1) + qJDD(3);
t88 = qJDD(4) + t91;
t99 = sin(qJ(4));
t110 = -t103 * t87 - t99 * t88;
t64 = t103 * t88 - t99 * t87;
t40 = -t100 * t64 + t104 * t110;
t41 = t100 * t110 + t104 * t64;
t102 = cos(qJ(5));
t101 = sin(qJ(1));
t105 = cos(qJ(1));
t112 = -t105 * g(2) - t101 * g(3);
t75 = qJDD(1) * pkin(1) + t112;
t107 = qJD(1) ^ 2;
t111 = t101 * g(2) - t105 * g(3);
t76 = -t107 * pkin(1) - t111;
t96 = sin(pkin(9));
t97 = cos(pkin(9));
t118 = t97 * t75 - t96 * t76;
t44 = qJDD(1) * pkin(2) + t118;
t124 = t96 * t75 + t97 * t76;
t45 = -t107 * pkin(2) + t124;
t29 = -t100 * t45 + t104 * t44;
t27 = t91 * pkin(3) + t29;
t30 = t100 * t44 + t104 * t45;
t90 = t92 ^ 2;
t28 = -t90 * pkin(3) + t30;
t21 = t103 * t28 + t99 * t27;
t19 = -t87 * pkin(4) + t88 * pkin(8) + t21;
t95 = -g(1) + qJDD(2);
t98 = sin(qJ(5));
t13 = -t102 * t95 + t98 * t19;
t14 = t102 * t19 + t98 * t95;
t7 = t102 * t14 + t98 * t13;
t20 = t103 * t27 - t99 * t28;
t18 = -t88 * pkin(4) - t87 * pkin(8) - t20;
t127 = -pkin(4) * t18 + pkin(8) * t7;
t79 = t102 * t87 * t98;
t73 = qJDD(5) + t79;
t126 = t98 * t73;
t125 = t98 * t88;
t74 = qJDD(5) - t79;
t123 = t102 * t74;
t122 = qJD(5) * t89;
t4 = -t103 * t18 + t99 * t7;
t121 = pkin(3) * t4 + t127;
t106 = qJD(5) ^ 2;
t93 = t98 ^ 2;
t82 = t93 * t87;
t77 = -t82 - t106;
t53 = -t98 * t77 - t123;
t59 = 0.2e1 * t102 * t122 + t125;
t120 = -pkin(4) * t59 + pkin(8) * t53 + t98 * t18;
t94 = t102 ^ 2;
t83 = t94 * t87;
t78 = -t83 - t106;
t51 = t102 * t78 - t126;
t81 = t102 * t88;
t60 = -0.2e1 * t98 * t122 + t81;
t119 = pkin(4) * t60 + pkin(8) * t51 - t102 * t18;
t62 = (t93 + t94) * t88;
t67 = t82 + t83;
t117 = pkin(4) * t67 + pkin(8) * t62 + t7;
t34 = -t103 * t59 + t99 * t53;
t116 = pkin(3) * t34 + t120;
t33 = t103 * t60 + t99 * t51;
t115 = pkin(3) * t33 + t119;
t114 = pkin(3) * t64 + t20;
t39 = t103 * t67 + t99 * t62;
t113 = pkin(3) * t39 + t117;
t71 = -t100 * t91 - t104 * t90;
t109 = t100 * t90 - t104 * t91;
t108 = pkin(3) * t110 - t21;
t52 = t123 + t98 * (t83 - t106);
t50 = t102 * (-t82 + t106) + t126;
t47 = t60 * t102;
t46 = t59 * t98;
t42 = t103 * t62 - t99 * t67;
t37 = t102 * t59 + t98 * t60;
t36 = t103 * t53 + t99 * t59;
t35 = t103 * t51 - t99 * t60;
t25 = t100 * t42 + t104 * t39;
t24 = t100 * t36 + t104 * t34;
t23 = t100 * t35 + t104 * t33;
t22 = t100 * t30 + t104 * t29;
t10 = t103 * t21 - t99 * t20;
t9 = t103 * t20 + t99 * t21;
t8 = pkin(3) * t9;
t5 = t103 * t7 + t99 * t18;
t2 = t100 * t10 + t104 * t9;
t1 = t100 * t5 + t104 * t4;
t3 = [0, 0, 0, 0, 0, qJDD(1), t112, t111, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t97 * qJDD(1) - t96 * t107) + t118, pkin(1) * (-t96 * qJDD(1) - t97 * t107) - t124, 0, pkin(1) * (t97 * t118 + t96 * t124), 0, 0, 0, 0, 0, t91, pkin(1) * (-t109 * t97 + t96 * t71) - pkin(2) * t109 + t29, pkin(1) * (t96 * t109 + t97 * t71) + pkin(2) * t71 - t30, 0, pkin(1) * (t96 * (-t100 * t29 + t104 * t30) + t97 * t22) + pkin(2) * t22, 0, 0, 0, 0, 0, t88, pkin(1) * (t96 * t40 + t97 * t41) + pkin(2) * t41 + t114, pkin(1) * (t97 * t40 - t96 * t41) + pkin(2) * t40 + t108, 0, pkin(1) * (t96 * (t104 * t10 - t100 * t9) + t97 * t2) + pkin(2) * t2 + t8, t46, t37, t50, t47, t52, 0, pkin(1) * (t96 * (-t100 * t33 + t104 * t35) + t97 * t23) + pkin(2) * t23 + t115, pkin(1) * (t96 * (-t100 * t34 + t104 * t36) + t97 * t24) + pkin(2) * t24 + t116, pkin(1) * (t96 * (-t100 * t39 + t104 * t42) + t97 * t25) + pkin(2) * t25 + t113, pkin(1) * (t96 * (-t100 * t4 + t104 * t5) + t97 * t1) + pkin(2) * t1 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, t102 * t73 + t98 * t78, t102 * t77 - t98 * t74, 0, -t102 * t13 + t98 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t29, -t30, 0, 0, 0, 0, 0, 0, 0, t88, t114, t108, 0, t8, t46, t37, t50, t47, t52, 0, t115, t116, t113, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t20, -t21, 0, 0, t46, t37, t50, t47, t52, 0, t119, t120, t117, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t82 - t83, t125, t79, t81, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t3;
