% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:47
% EndTime: 2019-12-05 15:28:49
% DurationCPUTime: 0.84s
% Computational Cost: add. (1272->197), mult. (2737->227), div. (0->0), fcn. (2010->8), ass. (0->108)
t139 = qJDD(2) * pkin(2);
t91 = pkin(7) + qJ(2);
t85 = sin(t91);
t87 = cos(t91);
t158 = g(1) * t85 - g(2) * t87;
t110 = t158 - qJDD(3) + t139;
t131 = qJD(2) * qJD(3);
t156 = qJ(3) * qJDD(2) + t131;
t119 = g(1) * t87 + g(2) * t85;
t149 = cos(qJ(4));
t122 = qJDD(2) * t149;
t93 = sin(pkin(8));
t134 = t93 * qJDD(2);
t94 = cos(pkin(8));
t96 = sin(qJ(4));
t117 = -t94 * t122 + t96 * t134;
t56 = t149 * t93 + t96 * t94;
t53 = t56 * qJD(4);
t24 = qJD(2) * t53 + t117;
t126 = t149 * t94;
t120 = qJD(2) * t126;
t145 = t96 * t93;
t127 = qJD(2) * t145;
t48 = -t120 + t127;
t123 = qJD(4) * t149;
t140 = qJD(4) * t96;
t52 = -t94 * t123 + t93 * t140;
t142 = -t56 * t24 + t52 * t48;
t107 = t126 - t145;
t157 = t56 * qJD(2);
t133 = t94 * qJDD(2);
t128 = qJD(4) * t120 + t93 * t122 + t96 * t133;
t23 = qJD(4) * t127 - t128;
t161 = t107 * t23 + t157 * t53;
t163 = t161 - t142;
t162 = t161 + t142;
t153 = t157 ^ 2;
t45 = t48 ^ 2;
t160 = -t45 - t153;
t159 = -t45 + t153;
t144 = pkin(6) + qJ(3);
t65 = t144 * t93;
t66 = t144 * t94;
t108 = -t149 * t65 - t96 * t66;
t17 = t107 * qJD(3) + t108 * qJD(4);
t35 = t149 * t66 - t96 * t65;
t90 = pkin(8) + qJ(4);
t84 = sin(t90);
t155 = -t17 * qJD(4) - t35 * qJDD(4) - t158 * t84;
t81 = t94 * qJDD(1);
t36 = t81 + (-t144 * qJDD(2) - t131) * t93;
t42 = t93 * qJDD(1) + t156 * t94;
t37 = pkin(6) * t133 + t42;
t83 = t94 * qJD(1);
t43 = -qJD(2) * t65 + t83;
t129 = t43 * t123 + t149 * t37 + t96 * t36;
t86 = cos(t90);
t154 = -g(3) * t84 - t119 * t86 + t129;
t148 = t157 * t48;
t135 = qJ(3) * qJD(2);
t58 = t93 * qJD(1) + t94 * t135;
t44 = t94 * qJD(2) * pkin(6) + t58;
t146 = t96 * t44;
t88 = t93 ^ 2;
t89 = t94 ^ 2;
t141 = t88 + t89;
t138 = qJDD(4) * pkin(4);
t20 = t149 * t44 + t96 * t43;
t137 = t20 * qJD(4);
t19 = t149 * t43 - t146;
t136 = qJD(5) - t19;
t130 = qJDD(4) * qJ(5);
t78 = t94 * pkin(3) + pkin(2);
t124 = t19 + t146;
t121 = t44 * t123 + t43 * t140 - t149 * t36 + t96 * t37;
t116 = t86 * pkin(4) + t84 * qJ(5);
t113 = -t107 * t24 + t48 * t53;
t112 = (-t93 * t135 + t83) * t93 - t58 * t94;
t27 = qJD(4) * t53 - qJDD(4) * t107;
t64 = -t78 * qJD(2) + qJD(3);
t60 = -t78 * qJDD(2) + qJDD(3);
t106 = -g(3) * t86 + t119 * t84 - t121;
t104 = t110 + t139;
t41 = -t156 * t93 + t81;
t103 = -t41 * t93 + t42 * t94 - t119;
t18 = t56 * qJD(3) + t35 * qJD(4);
t102 = -t18 * qJD(4) + t108 * qJDD(4) + t158 * t86;
t101 = t108 * t23 + t157 * t18 - t17 * t48 - t35 * t24 - t119;
t15 = t48 * pkin(4) - qJ(5) * t157 + t64;
t100 = t15 * t157 + qJDD(5) - t106;
t99 = t24 * pkin(4) + t23 * qJ(5) + t60;
t98 = 0.2e1 * t157 * qJD(4) + t117;
t92 = qJDD(1) - g(3);
t59 = t87 * t78;
t25 = -t52 * qJD(4) + t56 * qJDD(4);
t22 = -pkin(4) * t107 - t56 * qJ(5) - t78;
t21 = pkin(4) * t157 + t48 * qJ(5);
t16 = qJD(4) * qJ(5) + t20;
t14 = -qJD(4) * pkin(4) + t136;
t11 = t53 * pkin(4) + t52 * qJ(5) - t56 * qJD(5);
t10 = (t48 - t127) * qJD(4) + t128;
t9 = (t48 + t127) * qJD(4) - t128;
t6 = -t157 * t52 - t23 * t56;
t5 = -qJD(5) * t157 + t99;
t3 = -t44 * t140 + t129;
t2 = qJDD(5) + t121 - t138;
t1 = t130 + (qJD(5) - t146) * qJD(4) + t129;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t94 + t42 * t93 - g(3), 0, 0, 0, 0, 0, 0, -t27, -t25, t162, -t107 * t121 - t19 * t53 - t20 * t52 + t3 * t56 - g(3), 0, 0, 0, 0, 0, 0, -t27, t162, t25, t1 * t56 - t107 * t2 + t14 * t53 - t16 * t52 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t158, t119, 0, 0, t88 * qJDD(2), 0.2e1 * t93 * t133, 0, t89 * qJDD(2), 0, 0, t104 * t94, -t104 * t93, t156 * t141 + t103, t110 * pkin(2) + t103 * qJ(3) - t112 * qJD(3), t6, -t163, t25, t113, -t27, 0, -t107 * t60 - t78 * t24 + t64 * t53 + t102, t78 * t23 - t64 * t52 + t60 * t56 + t155, t107 * t3 + t121 * t56 + t19 * t52 - t20 * t53 + t101, t3 * t35 + t20 * t17 - t121 * t108 - t19 * t18 - t60 * t78 - g(1) * (t144 * t87 - t85 * t78) - g(2) * (t144 * t85 + t59), t6, t25, t163, 0, t27, t113, -t107 * t5 + t11 * t48 + t15 * t53 + t22 * t24 + t102, t1 * t107 - t14 * t52 - t16 * t53 + t2 * t56 + t101, -t11 * t157 + t15 * t52 + t22 * t23 - t5 * t56 - t155, -g(2) * t59 + t1 * t35 + t15 * t11 + t14 * t18 + t16 * t17 - t2 * t108 + t5 * t22 + (-g(1) * t144 - g(2) * t116) * t87 + (-g(1) * (-t116 - t78) - g(2) * t144) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t134, -t141 * qJD(2) ^ 2, t112 * qJD(2) - t110, 0, 0, 0, 0, 0, 0, t98, -t9, t160, t157 * t19 + t20 * t48 - t158 + t60, 0, 0, 0, 0, 0, 0, t98, t160, t9, t16 * t48 + (-qJD(5) - t14) * t157 + t99 - t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t159, t10, -t148, -t117, qJDD(4), -t157 * t64 + t106 + t137, t124 * qJD(4) + t64 * t48 - t154, 0, 0, t148, t10, -t159, qJDD(4), t117, -t148, -t21 * t48 - t100 + t137 + 0.2e1 * t138, pkin(4) * t23 - t24 * qJ(5) + (t16 - t20) * t157 + (t14 - t136) * t48, 0.2e1 * t130 - t15 * t48 + t21 * t157 + (0.2e1 * qJD(5) - t124) * qJD(4) + t154, -t2 * pkin(4) - g(3) * t116 + t1 * qJ(5) + t136 * t16 - t14 * t20 - t15 * t21 + t119 * (pkin(4) * t84 - qJ(5) * t86); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t148, t10, -qJD(4) ^ 2 - t153, -t16 * qJD(4) + t100 - t138;];
tau_reg = t4;
