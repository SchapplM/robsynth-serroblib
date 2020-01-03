% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:52
% DurationCPUTime: 1.02s
% Computational Cost: add. (3078->195), mult. (5563->258), div. (0->0), fcn. (2868->8), ass. (0->126)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t103 = sin(qJ(4));
t131 = qJD(1) * t103;
t66 = t105 * qJD(4) + t102 * t131;
t67 = -t102 * qJD(4) + t105 * t131;
t51 = t66 * t67;
t106 = cos(qJ(4));
t128 = t106 * qJDD(1);
t130 = qJD(1) * qJD(4);
t85 = t103 * t130;
t72 = t85 - t128;
t65 = qJDD(5) - t72;
t146 = -t51 + t65;
t148 = t102 * t146;
t147 = t105 * t146;
t123 = t106 * t130;
t129 = t103 * qJDD(1);
t70 = -t123 - t129;
t122 = -t105 * qJDD(4) + t102 * t70;
t83 = t106 * qJD(1) + qJD(5);
t32 = (qJD(5) - t83) * t67 - t122;
t63 = t66 ^ 2;
t64 = t67 ^ 2;
t81 = t83 ^ 2;
t109 = qJD(1) ^ 2;
t145 = pkin(1) + pkin(2);
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t144 = t104 * g(1) - t107 * g(2);
t126 = -qJDD(2) + t144;
t114 = -t109 * qJ(2) - t126;
t110 = -t145 * qJDD(1) + t114;
t120 = t107 * g(1) + t104 * g(2);
t115 = 0.2e1 * qJD(2) * qJD(1) - t120;
t94 = qJDD(1) * qJ(2);
t112 = t115 + t94;
t57 = -t145 * t109 + t112;
t98 = sin(pkin(8));
t99 = cos(pkin(8));
t38 = t98 * t110 + t99 * t57;
t108 = qJD(4) ^ 2;
t121 = pkin(4) * t106 + pkin(7) * t103;
t113 = t109 * t121;
t134 = g(3) + qJDD(3);
t30 = -t109 * pkin(3) - qJDD(1) * pkin(6) + t38;
t27 = t103 * t30 - t106 * t134;
t21 = -qJDD(4) * pkin(4) - t108 * pkin(7) - t103 * t113 + t27;
t143 = t102 * t21;
t43 = t51 + t65;
t142 = t102 * t43;
t141 = t102 * t83;
t82 = t106 * t109 * t103;
t77 = qJDD(4) + t82;
t140 = t103 * t77;
t139 = t105 * t21;
t138 = t105 * t43;
t137 = t105 * t83;
t78 = qJDD(4) - t82;
t136 = t106 * t78;
t135 = qJDD(1) * pkin(1);
t132 = qJD(5) + t83;
t127 = t106 * t51;
t125 = qJ(2) * t99 - pkin(6);
t118 = -t70 + t123;
t119 = -t72 - t85;
t124 = -t99 * t110 + t98 * t57;
t29 = qJDD(1) * pkin(3) - t109 * pkin(6) + t124;
t20 = t119 * pkin(4) + t118 * pkin(7) + t29;
t28 = t103 * t134 + t106 * t30;
t22 = -t108 * pkin(4) + qJDD(4) * pkin(7) - t106 * t113 + t28;
t10 = t102 * t22 - t105 * t20;
t11 = t102 * t20 + t105 * t22;
t4 = t102 * t10 + t105 * t11;
t117 = t105 * t10 - t102 * t11;
t14 = t103 * t27 + t106 * t28;
t116 = -t102 * qJDD(4) - t105 * t70;
t111 = qJ(2) * t98 + pkin(3) + t121;
t46 = t66 * qJD(5) - t116;
t96 = t106 ^ 2;
t95 = t103 ^ 2;
t90 = t96 * t109;
t89 = t95 * t109;
t80 = -t90 - t108;
t79 = -t89 - t108;
t76 = t89 + t90;
t75 = (-t95 - t96) * qJDD(1);
t74 = t99 * qJDD(1) + t98 * t109;
t73 = -t98 * qJDD(1) + t99 * t109;
t71 = 0.2e1 * t85 - t128;
t69 = 0.2e1 * t123 + t129;
t61 = -t114 + t135;
t60 = t83 * t66;
t59 = -t64 + t81;
t58 = t63 - t81;
t54 = -t103 * t79 - t136;
t53 = t106 * t80 - t140;
t50 = t64 - t63;
t49 = -t64 - t81;
t48 = t98 * t75 + t99 * t76;
t47 = -t81 - t63;
t45 = t67 * qJD(5) - t122;
t41 = t63 + t64;
t40 = t98 * t54 + t99 * t69;
t39 = t98 * t53 + t99 * t71;
t36 = -t132 * t66 + t116;
t35 = t46 - t60;
t34 = t46 + t60;
t33 = t132 * t67 - t122;
t26 = -t102 * t49 - t138;
t25 = t105 * t49 - t142;
t24 = t105 * t47 - t148;
t23 = t102 * t47 + t147;
t18 = -t124 * t99 + t98 * t38;
t17 = t102 * t35 + t105 * t32;
t16 = t102 * t32 - t105 * t35;
t15 = -t103 * t36 + t106 * t26;
t13 = -t103 * t33 + t106 * t24;
t12 = -t103 * t41 + t106 * t17;
t8 = t98 * t14 - t99 * t29;
t7 = t98 * t15 - t99 * t25;
t6 = t98 * t13 - t99 * t23;
t5 = t98 * t12 - t99 * t16;
t2 = t103 * t21 + t106 * t4;
t1 = t117 * t99 + t98 * t2;
t3 = [0, 0, 0, 0, 0, qJDD(1), t144, t120, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t126 + 0.2e1 * t135, 0, t115 + 0.2e1 * t94, pkin(1) * t61 + qJ(2) * (-t109 * pkin(1) + t112), 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t73 + t145 * t74 + t124, qJ(2) * t74 + t145 * t73 + t38, 0, qJ(2) * (t124 * t98 + t99 * t38) - t145 * t18, t118 * t103, -t103 * t71 + t106 * t69, -t140 - t106 * (-t89 + t108), t119 * t106, -t103 * (t90 - t108) - t136, 0, qJ(2) * (t99 * t53 - t98 * t71) + t106 * t29 - pkin(3) * t71 - pkin(6) * t53 - t145 * t39, qJ(2) * (t99 * t54 - t98 * t69) - t103 * t29 - pkin(3) * t69 - pkin(6) * t54 - t145 * t40, qJ(2) * (t99 * t75 - t98 * t76) - pkin(3) * t76 - pkin(6) * t75 - t145 * t48 - t14, qJ(2) * (t99 * t14 + t98 * t29) + pkin(3) * t29 - pkin(6) * t14 - t145 * t8, -t103 * (t105 * t46 + t67 * t141) + t127, -t103 * (-t102 * t34 + t105 * t33) + t106 * t50, -t103 * (-t102 * t59 + t147) + t106 * t35, -t103 * (-t102 * t45 - t66 * t137) - t127, -t103 * (t105 * t58 - t142) + t106 * t32, t106 * t65 - t103 * (-t102 * t67 + t105 * t66) * t83, qJ(2) * (t99 * t13 + t98 * t23) - t103 * (-pkin(7) * t23 + t143) - t106 * (-pkin(4) * t23 + t10) + pkin(3) * t23 - pkin(6) * t13 - t145 * t6, qJ(2) * (t99 * t15 + t98 * t25) - t103 * (-pkin(7) * t25 + t139) - t106 * (-pkin(4) * t25 + t11) + pkin(3) * t25 - pkin(6) * t15 - t145 * t7, -t103 * t117 + t111 * t16 + t125 * t12 - t145 * t5, -t145 * t1 - t111 * t117 + t125 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t109, -t61, 0, 0, 0, 0, 0, 0, -t74, -t73, 0, t18, 0, 0, 0, 0, 0, 0, t39, t40, t48, t8, 0, 0, 0, 0, 0, 0, t6, t7, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, 0, 0, 0, 0, t103 * t80 + t106 * t77, -t103 * t78 + t106 * t79, 0, t103 * t28 - t106 * t27, 0, 0, 0, 0, 0, 0, t103 * t24 + t106 * t33, t103 * t26 + t106 * t36, t103 * t17 + t106 * t41, t103 * t4 - t106 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t89 - t90, -t129, t82, -t128, qJDD(4), -t27, -t28, 0, 0, t102 * t46 - t67 * t137, t102 * t33 + t105 * t34, t105 * t59 + t148, t105 * t45 - t66 * t141, t102 * t58 + t138, (t102 * t66 + t105 * t67) * t83, pkin(4) * t33 + pkin(7) * t24 - t139, pkin(4) * t36 + pkin(7) * t26 + t143, pkin(4) * t41 + pkin(7) * t17 + t4, -pkin(4) * t21 + pkin(7) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t50, t35, -t51, t32, t65, -t10, -t11, 0, 0;];
tauJ_reg = t3;
