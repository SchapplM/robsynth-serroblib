% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:06
% EndTime: 2019-12-31 17:44:07
% DurationCPUTime: 0.82s
% Computational Cost: add. (984->189), mult. (1898->239), div. (0->0), fcn. (1316->10), ass. (0->114)
t89 = sin(pkin(8));
t86 = t89 ^ 2;
t91 = cos(pkin(8));
t87 = t91 ^ 2;
t151 = t86 + t87;
t124 = qJD(1) * qJD(4);
t117 = t89 * t124;
t133 = t89 * qJ(4);
t116 = pkin(2) + t133;
t92 = cos(pkin(7));
t145 = t92 * pkin(1);
t36 = (pkin(3) + pkin(4)) * t91 + t116 + t145;
t10 = t36 * qJDD(1) - qJDD(3) + t117;
t88 = qJ(1) + pkin(7);
t82 = sin(t88);
t83 = cos(t88);
t134 = g(1) * t83 + g(2) * t82;
t132 = pkin(1) * qJDD(1);
t93 = sin(qJ(5));
t95 = cos(qJ(5));
t52 = t89 * t93 + t91 * t95;
t43 = t52 * qJD(1);
t53 = t89 * t95 - t91 * t93;
t130 = qJD(1) * t91;
t122 = t93 * t130;
t101 = qJD(5) * t122 - t52 * qJDD(1);
t149 = t43 ^ 2;
t131 = qJD(1) * t89;
t121 = t95 * t131;
t45 = t121 - t122;
t148 = t45 ^ 2;
t146 = g(1) * t82;
t73 = g(2) * t83;
t90 = sin(pkin(7));
t68 = t90 * pkin(1) + qJ(3);
t144 = -pkin(6) + t68;
t12 = qJD(5) * t121 - t101;
t41 = t52 * qJD(5);
t143 = -t53 * t12 + t41 * t43;
t125 = qJD(1) * qJD(3);
t49 = qJDD(1) * t68 + t125;
t29 = t89 * qJDD(2) + t91 * t49;
t142 = t29 * t89 - g(3);
t23 = t29 * t91;
t58 = t68 * qJD(1);
t38 = t89 * qJD(2) + t91 * t58;
t141 = t38 * t91;
t140 = t45 * t43;
t139 = t82 * t91;
t138 = t83 * t91;
t135 = qJD(3) * t141 + t68 * t23;
t129 = qJD(4) * t89;
t105 = -t91 * pkin(3) - t116;
t46 = t105 - t145;
t128 = qJDD(1) * t46;
t76 = -pkin(2) - t145;
t127 = qJDD(1) * t76;
t77 = t86 * qJDD(1);
t80 = t87 * qJDD(1);
t79 = t89 * qJDD(1);
t126 = t91 * qJDD(1);
t96 = cos(qJ(1));
t123 = t96 * pkin(1) + t83 * pkin(2) + t82 * qJ(3);
t120 = t89 * t126;
t94 = sin(qJ(1));
t119 = -t94 * pkin(1) + t83 * qJ(3);
t118 = t73 - t146;
t28 = t91 * qJDD(2) - t89 * t49;
t26 = qJDD(4) - t28;
t18 = -pkin(6) * t79 + t26;
t20 = -pkin(6) * t126 + t29;
t115 = t95 * t18 - t93 * t20;
t37 = t91 * qJD(2) - t89 * t58;
t114 = pkin(3) * t138 + t83 * t133 + t123;
t112 = g(1) * t94 - g(2) * t96;
t111 = t93 * t126 - t95 * t79;
t11 = qJD(1) * t41 + t111;
t42 = t53 * qJD(5);
t110 = -t52 * t11 + t45 * t42;
t109 = t93 * t18 + t95 * t20;
t35 = qJD(4) - t37;
t24 = -pkin(6) * t131 + t35;
t25 = -pkin(6) * t130 + t38;
t3 = t95 * t24 - t93 * t25;
t4 = t93 * t24 + t95 * t25;
t47 = t144 * t89;
t48 = t144 * t91;
t8 = t95 * t47 - t93 * t48;
t9 = t93 * t47 + t95 * t48;
t107 = -t82 * pkin(2) + t119;
t100 = qJDD(3) + t128;
t19 = t100 - t117;
t104 = -t128 - t19 - t73;
t56 = qJDD(3) + t127;
t103 = -t127 - t56 - t73;
t102 = -t134 + t23 + (t77 + t80) * t68 + t151 * t125;
t98 = qJD(1) ^ 2;
t97 = qJD(5) ^ 2;
t84 = g(3) * t91;
t61 = g(1) * t139;
t57 = t151 * t98;
t34 = t52 * t83;
t33 = t53 * t83;
t32 = t52 * t82;
t31 = t53 * t82;
t30 = t46 * qJD(1) + qJD(3);
t21 = t36 * qJD(1) - qJD(3);
t14 = -t42 * qJD(5) - t52 * qJDD(5);
t13 = -t41 * qJD(5) + t53 * qJDD(5);
t6 = t53 * qJD(3) - t9 * qJD(5);
t5 = t52 * qJD(3) + t8 * qJD(5);
t2 = -t4 * qJD(5) + t115;
t1 = t3 * qJD(5) + t109;
t7 = [0, 0, 0, 0, 0, qJDD(1), t112, g(1) * t96 + g(2) * t94, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t92 * t132 - t118, -0.2e1 * t90 * t132 + t134, 0, (t112 + (t90 ^ 2 + t92 ^ 2) * t132) * pkin(1), t77, 0.2e1 * t120, 0, t80, 0, 0, t103 * t91 + t61, (-t103 - t146) * t89, -t28 * t89 + t102, t56 * t76 - g(1) * t107 - g(2) * t123 + (-t37 * qJD(3) - t28 * t68) * t89 + t135, t77, 0, -0.2e1 * t120, 0, 0, t80, t61 + (t104 + t117) * t91, t26 * t89 + t102, t86 * t124 + (t104 + t146) * t89, t19 * t46 - g(1) * (-pkin(3) * t139 + t107) - g(2) * t114 + (qJ(4) * t146 + t35 * qJD(3) - t30 * qJD(4) + t26 * t68) * t89 + t135, -t11 * t53 - t45 * t41, -t110 + t143, t13, t12 * t52 + t43 * t42, t14, 0, g(1) * t32 - g(2) * t34 + t6 * qJD(5) + t8 * qJDD(5) + t10 * t52 + t36 * t12 + t43 * t129 + t21 * t42, g(1) * t31 - g(2) * t33 - t5 * qJD(5) - t9 * qJDD(5) + t10 * t53 - t36 * t11 + t45 * t129 - t21 * t41, -t1 * t52 + t8 * t11 - t9 * t12 - t2 * t53 + t3 * t41 - t4 * t42 - t5 * t43 - t6 * t45 + t134, t1 * t9 + t4 * t5 + t2 * t8 + t3 * t6 + t10 * t36 + t21 * t129 - g(1) * (-t83 * pkin(6) + t119) - g(2) * (pkin(4) * t138 + t114) + (-g(1) * (-t91 * pkin(4) + t105) + g(2) * pkin(6)) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t91 + t142, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t91 + t142, 0, 0, 0, 0, 0, 0, t14, -t13, t110 + t143, t1 * t53 - t2 * t52 - t3 * t42 - t4 * t41 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t79, -t57, (t37 * t89 - t141) * qJD(1) + t56 + t118, 0, 0, 0, 0, 0, 0, -t126, -t57, -t79, (-t141 + (-qJD(4) - t35) * t89) * qJD(1) + t100 + t118, 0, 0, 0, 0, 0, 0, (-t45 - t121) * qJD(5) + t101, 0.2e1 * t43 * qJD(5) + t111, t148 + t149, -t3 * t45 - t4 * t43 - t10 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * t98 * t91, t79, -t86 * t98, t84 + (qJD(1) * t30 - t134) * t89 + t26, 0, 0, 0, 0, 0, 0, qJDD(5) * t95 - t43 * t131 - t97 * t93, -qJDD(5) * t93 - t45 * t131 - t97 * t95, t95 * t11 - t93 * t12 + (-t43 * t95 + t45 * t93) * qJD(5), t1 * t93 + t2 * t95 + t84 + (-t3 * t93 + t4 * t95) * qJD(5) + (-qJD(1) * t21 - t134) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t148 - t149, -t111, -t140, (t45 - t121) * qJD(5) + t101, qJDD(5), -g(1) * t33 - g(2) * t31 + g(3) * t52 - t21 * t45 + t115, g(1) * t34 + g(2) * t32 + g(3) * t53 + t21 * t43 - t109, 0, 0;];
tau_reg = t7;
