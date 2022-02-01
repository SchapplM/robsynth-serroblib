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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:34:45
% EndTime: 2022-01-23 09:34:47
% DurationCPUTime: 0.71s
% Computational Cost: add. (3491->125), mult. (4975->187), div. (0->0), fcn. (2890->10), ass. (0->97)
t101 = sin(qJ(3));
t105 = cos(qJ(3));
t100 = sin(qJ(4));
t104 = cos(qJ(4));
t93 = qJD(1) + qJD(3);
t89 = qJD(4) + t93;
t87 = t89 ^ 2;
t92 = qJDD(1) + qJDD(3);
t88 = qJDD(4) + t92;
t111 = -t100 * t88 - t104 * t87;
t65 = t100 * t87 - t104 * t88;
t40 = t101 * t65 + t105 * t111;
t131 = -t101 * t111 + t105 * t65;
t103 = cos(qJ(5));
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t119 = t102 * g(1) - t106 * g(2);
t75 = qJDD(1) * pkin(1) + t119;
t108 = qJD(1) ^ 2;
t112 = t106 * g(1) + t102 * g(2);
t76 = -t108 * pkin(1) - t112;
t97 = sin(pkin(9));
t98 = cos(pkin(9));
t118 = t98 * t75 - t97 * t76;
t44 = qJDD(1) * pkin(2) + t118;
t125 = t97 * t75 + t98 * t76;
t45 = -t108 * pkin(2) + t125;
t29 = -t101 * t45 + t105 * t44;
t27 = t92 * pkin(3) + t29;
t30 = t101 * t44 + t105 * t45;
t91 = t93 ^ 2;
t28 = -t91 * pkin(3) + t30;
t21 = t100 * t27 + t104 * t28;
t19 = -t87 * pkin(4) + t88 * pkin(8) + t21;
t96 = -g(3) + qJDD(2);
t99 = sin(qJ(5));
t13 = -t103 * t96 + t99 * t19;
t14 = t103 * t19 + t99 * t96;
t7 = t103 * t14 + t99 * t13;
t20 = -t100 * t28 + t104 * t27;
t18 = -t88 * pkin(4) - t87 * pkin(8) - t20;
t128 = -pkin(4) * t18 + pkin(8) * t7;
t79 = t103 * t87 * t99;
t73 = qJDD(5) + t79;
t127 = t99 * t73;
t126 = t99 * t88;
t74 = qJDD(5) - t79;
t124 = t103 * t74;
t123 = qJD(5) * t89;
t4 = t100 * t7 - t104 * t18;
t122 = pkin(3) * t4 + t128;
t107 = qJD(5) ^ 2;
t94 = t99 ^ 2;
t82 = t94 * t87;
t77 = -t82 - t107;
t53 = -t99 * t77 - t124;
t59 = 0.2e1 * t103 * t123 + t126;
t121 = -pkin(4) * t59 + pkin(8) * t53 + t99 * t18;
t95 = t103 ^ 2;
t83 = t95 * t87;
t78 = -t83 - t107;
t52 = t103 * t78 - t127;
t81 = t103 * t88;
t60 = -0.2e1 * t123 * t99 + t81;
t120 = pkin(4) * t60 + pkin(8) * t52 - t103 * t18;
t62 = (t94 + t95) * t88;
t67 = t82 + t83;
t117 = pkin(4) * t67 + pkin(8) * t62 + t7;
t34 = t100 * t53 - t104 * t59;
t116 = pkin(3) * t34 + t121;
t33 = t100 * t52 + t104 * t60;
t115 = pkin(3) * t33 + t120;
t114 = -pkin(3) * t65 + t20;
t39 = t100 * t62 + t104 * t67;
t113 = pkin(3) * t39 + t117;
t71 = -t101 * t92 - t105 * t91;
t110 = t101 * t91 - t105 * t92;
t109 = pkin(3) * t111 - t21;
t51 = t127 + t103 * (-t82 + t107);
t50 = t99 * (t83 - t107) + t124;
t47 = t59 * t99;
t46 = t60 * t103;
t42 = -t100 * t67 + t104 * t62;
t37 = t103 * t59 + t99 * t60;
t36 = t100 * t59 + t104 * t53;
t35 = -t100 * t60 + t104 * t52;
t25 = t101 * t42 + t105 * t39;
t24 = t101 * t36 + t105 * t34;
t23 = t101 * t35 + t105 * t33;
t22 = t101 * t30 + t105 * t29;
t10 = -t100 * t20 + t104 * t21;
t9 = t100 * t21 + t104 * t20;
t8 = pkin(3) * t9;
t5 = t100 * t18 + t104 * t7;
t2 = t101 * t10 + t105 * t9;
t1 = t101 * t5 + t105 * t4;
t3 = [0, 0, 0, 0, 0, qJDD(1), t119, t112, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t98 * qJDD(1) - t97 * t108) + t118, pkin(1) * (-t97 * qJDD(1) - t98 * t108) - t125, 0, pkin(1) * (t98 * t118 + t97 * t125), 0, 0, 0, 0, 0, t92, pkin(1) * (-t110 * t98 + t71 * t97) - pkin(2) * t110 + t29, pkin(1) * (t110 * t97 + t98 * t71) + pkin(2) * t71 - t30, 0, pkin(1) * (t97 * (-t101 * t29 + t105 * t30) + t98 * t22) + pkin(2) * t22, 0, 0, 0, 0, 0, t88, pkin(1) * (-t131 * t98 + t97 * t40) - pkin(2) * t131 + t114, pkin(1) * (t97 * t131 + t98 * t40) + pkin(2) * t40 + t109, 0, pkin(1) * (t97 * (t105 * t10 - t101 * t9) + t98 * t2) + pkin(2) * t2 + t8, t47, t37, t51, t46, t50, 0, pkin(1) * (t97 * (-t101 * t33 + t105 * t35) + t98 * t23) + pkin(2) * t23 + t115, pkin(1) * (t97 * (-t101 * t34 + t105 * t36) + t98 * t24) + pkin(2) * t24 + t116, pkin(1) * (t97 * (-t101 * t39 + t105 * t42) + t98 * t25) + pkin(2) * t25 + t113, pkin(1) * (t97 * (-t101 * t4 + t105 * t5) + t98 * t1) + pkin(2) * t1 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, t103 * t73 + t99 * t78, t103 * t77 - t99 * t74, 0, -t103 * t13 + t99 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t29, -t30, 0, 0, 0, 0, 0, 0, 0, t88, t114, t109, 0, t8, t47, t37, t51, t46, t50, 0, t115, t116, t113, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t20, -t21, 0, 0, t47, t37, t51, t46, t50, 0, t120, t121, t117, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t82 - t83, t126, t79, t81, qJDD(5), -t13, -t14, 0, 0;];
tauJ_reg = t3;
