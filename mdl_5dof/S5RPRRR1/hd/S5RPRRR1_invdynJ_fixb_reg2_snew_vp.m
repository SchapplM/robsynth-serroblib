% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:51
% EndTime: 2019-12-05 18:09:58
% DurationCPUTime: 1.19s
% Computational Cost: add. (1808->197), mult. (3402->281), div. (0->0), fcn. (2668->8), ass. (0->135)
t103 = cos(qJ(5));
t100 = sin(qJ(4));
t104 = cos(qJ(4));
t101 = sin(qJ(3));
t130 = qJD(1) * t101;
t80 = t100 * qJD(3) + t104 * t130;
t105 = cos(qJ(3));
t89 = t105 * qJD(1) - qJD(4);
t99 = sin(qJ(5));
t62 = t103 * t89 + t99 * t80;
t64 = t103 * t80 - t99 * t89;
t41 = t64 * t62;
t129 = qJD(1) * qJD(3);
t121 = t105 * t129;
t92 = t101 * qJDD(1);
t82 = t92 + t121;
t119 = -t104 * qJDD(3) + t100 * t82;
t53 = -t80 * qJD(4) - t119;
t52 = qJDD(5) - t53;
t158 = -t41 + t52;
t160 = t158 * t99;
t159 = t103 * t158;
t108 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t106 = cos(qJ(1));
t123 = t102 * g(1) - t106 * g(2);
t116 = -qJDD(2) + t123;
t110 = -t108 * qJ(2) - t116;
t115 = t106 * g(1) + t102 * g(2);
t111 = 0.2e1 * qJD(2) * qJD(1) - t115;
t128 = qJDD(1) * qJ(2);
t72 = t111 + t128;
t66 = -t101 * g(3) + t105 * t72;
t44 = t100 * t66 - t104 * t110;
t42 = t104 * t44;
t45 = t100 * t110 + t104 * t66;
t157 = -t100 * t45 + t42;
t78 = -t104 * qJD(3) + t100 * t130;
t153 = t80 * t78;
t127 = t105 * qJDD(1);
t91 = t101 * t129;
t117 = -t91 + t127;
t77 = -qJDD(4) + t117;
t109 = -t77 - t153;
t156 = t100 * t109;
t155 = t104 * t109;
t36 = (qJD(4) + t89) * t80 + t119;
t112 = -t100 * qJDD(3) - t104 * t82;
t54 = -t78 * qJD(4) - t112;
t122 = t103 * t77 + t99 * t54;
t74 = qJD(5) + t78;
t17 = (qJD(5) - t74) * t64 + t122;
t60 = t62 ^ 2;
t61 = t64 ^ 2;
t73 = t74 ^ 2;
t75 = t78 ^ 2;
t76 = t80 ^ 2;
t87 = t89 ^ 2;
t154 = t74 * t99;
t24 = t41 + t52;
t152 = t99 * t24;
t97 = t101 ^ 2;
t98 = t105 ^ 2;
t151 = t97 + t98;
t150 = t100 * t44;
t49 = t77 - t153;
t148 = t100 * t49;
t147 = t100 * t89;
t65 = t105 * g(3) + t101 * t72;
t146 = t101 * t65;
t88 = t105 * t108 * t101;
t145 = t101 * (qJDD(3) + t88);
t144 = t103 * t24;
t143 = t103 * t44;
t142 = t103 * t74;
t141 = t104 * t49;
t140 = t104 * t89;
t139 = t105 * t44;
t138 = t105 * (qJDD(3) - t88);
t137 = t97 * t108;
t136 = t98 * t108;
t135 = t101 * t104;
t134 = qJD(4) - t89;
t131 = qJD(5) + t74;
t126 = t100 * t41;
t125 = t104 * t41;
t124 = t105 * t153;
t120 = t104 * t45 + t150;
t118 = t105 * t66 + t146;
t31 = -t103 * t65 + t99 * t45;
t32 = t103 * t45 + t99 * t65;
t114 = t103 * t32 + t99 * t31;
t6 = t103 * t31 - t99 * t32;
t113 = -t103 * t54 + t99 * t77;
t29 = -t62 * qJD(5) - t113;
t107 = qJD(3) ^ 2;
t83 = -0.2e1 * t91 + t127;
t81 = t92 + 0.2e1 * t121;
t69 = t78 * t89;
t68 = -t76 + t87;
t67 = t75 - t87;
t57 = t76 - t75;
t56 = -t76 - t87;
t55 = -t87 - t75;
t48 = t74 * t62;
t47 = -t61 + t73;
t46 = t60 - t73;
t40 = t54 - t69;
t39 = t54 + t69;
t37 = -t134 * t80 - t119;
t35 = t61 - t60;
t34 = -t61 - t73;
t33 = -t73 - t60;
t30 = t60 + t61;
t28 = -t64 * qJD(5) - t122;
t27 = (-t103 * t62 + t64 * t99) * t74;
t26 = (-t103 * t64 - t62 * t99) * t74;
t22 = t131 * t62 + t113;
t21 = t29 + t48;
t20 = t29 - t48;
t18 = -t131 * t64 - t122;
t16 = t103 * t29 - t64 * t154;
t15 = t64 * t142 + t99 * t29;
t14 = t62 * t142 - t99 * t28;
t13 = -t103 * t28 - t62 * t154;
t12 = t103 * t46 - t152;
t11 = -t99 * t47 + t159;
t10 = t99 * t46 + t144;
t9 = t103 * t47 + t160;
t8 = -t99 * t34 - t144;
t4 = t103 * t33 - t160;
t3 = -t103 * t17 + t99 * t21;
t2 = t103 * t18 - t99 * t20;
t1 = t103 * t20 + t99 * t18;
t5 = [0, 0, 0, 0, 0, qJDD(1), t123, t115, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t116, 0, t111 + 0.2e1 * t128, qJ(2) * t72, (t82 + t121) * t101, t101 * t83 + t105 * t81, t145 + t105 * (t107 - t137), (t117 - t91) * t105, t101 * (-t107 + t136) + t138, 0, -t105 * t110 + qJ(2) * (t105 * (-t107 - t136) - t145), t101 * t110 + qJ(2) * (-t138 - t101 * (-t107 - t137)), t151 * t128 + t118, qJ(2) * t118, t101 * (t104 * t54 + t80 * t147) - t124, t101 * (-t100 * t39 + t104 * t37) - t105 * t57, t101 * (-t100 * t68 + t155) - t105 * t40, t101 * (-t100 * t53 - t78 * t140) + t124, t101 * (t104 * t67 + t148) + t105 * t36, t105 * t77 + t101 * (-t100 * t80 + t104 * t78) * t89, t100 * t146 + t139 + qJ(2) * (t105 * (t104 * t55 - t156) - t101 * t37), t65 * t135 + t105 * t45 + qJ(2) * (t105 * (-t100 * t56 + t141) - t101 * (t134 * t78 + t112)), t101 * t157 + qJ(2) * (t105 * (t100 * t40 - t104 * t36) - t101 * (t75 + t76)), qJ(2) * (t105 * t120 + t146), t101 * (t104 * t16 + t126) - t105 * t15, t101 * (t100 * t35 + t104 * t2) - t105 * t1, t101 * (t100 * t21 + t104 * t11) - t105 * t9, t101 * (t104 * t14 - t126) + t105 * t13, t101 * (-t100 * t17 + t104 * t12) - t105 * t10, t101 * (t100 * t52 + t104 * t27) - t105 * t26, t101 * (-t100 * t31 + t99 * t42) + t103 * t139 + qJ(2) * (t105 * (-t100 * t18 + t104 * t4) - t101 * (-t99 * t33 - t159)), t101 * (-t100 * t32 + t103 * t42) - t99 * t139 + qJ(2) * (t105 * (-t100 * t22 + t104 * t8) - t101 * (-t103 * t34 + t152)), t6 * t135 - t105 * t114 + qJ(2) * (t105 * (-t100 * t30 + t104 * t3) - t101 * (t103 * t21 + t17 * t99)), qJ(2) * (t105 * (t104 * t114 + t150) - t101 * t6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t108, t110, 0, 0, 0, 0, 0, 0, -t83, t81, -t151 * t108, t110, 0, 0, 0, 0, 0, 0, t100 * t55 + t155, t104 * t56 + t148, -t100 * t36 - t104 * t40, -t157, 0, 0, 0, 0, 0, 0, t100 * t4 + t104 * t18, t100 * t8 + t104 * t22, t100 * t3 + t104 * t30, t100 * t114 - t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, (t97 - t98) * t108, t92, t88, t127, qJDD(3), -t65, -t66, 0, 0, t100 * t54 - t80 * t140, t100 * t37 + t104 * t39, t104 * t68 + t156, t104 * t53 - t78 * t147, t100 * t67 - t141, (t100 * t78 + t104 * t80) * t89, -t104 * t65, t100 * t65, t120, 0, t100 * t16 - t125, t100 * t2 - t104 * t35, t100 * t11 - t104 * t21, t100 * t14 + t125, t100 * t12 + t104 * t17, t100 * t27 - t104 * t52, t104 * t31 + t99 * t150, t100 * t143 + t104 * t32, t100 * t6, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t57, t40, -t153, -t36, -t77, -t44, -t45, 0, 0, t15, t1, t9, -t13, t10, t26, -t143, t99 * t44, t114, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t35, t21, -t41, -t17, t52, -t31, -t32, 0, 0;];
tauJ_reg = t5;
