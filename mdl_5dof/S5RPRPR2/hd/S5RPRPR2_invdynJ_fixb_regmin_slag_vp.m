% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:41
% EndTime: 2019-12-05 17:49:44
% DurationCPUTime: 0.68s
% Computational Cost: add. (1050->152), mult. (1737->196), div. (0->0), fcn. (1138->16), ass. (0->107)
t92 = qJ(1) + pkin(8);
t84 = qJ(3) + t92;
t74 = sin(t84);
t75 = cos(t84);
t158 = g(2) * t75 + g(3) * t74;
t94 = sin(pkin(8));
t154 = pkin(1) * t94;
t134 = qJD(3) * t154;
t96 = cos(pkin(8));
t76 = t96 * pkin(1) + pkin(2);
t160 = -qJD(1) * t134 + t76 * qJDD(1);
t62 = t76 * qJD(1);
t159 = qJD(3) * t62 + qJDD(1) * t154;
t100 = cos(qJ(5));
t95 = cos(pkin(9));
t140 = t100 * t95;
t93 = sin(pkin(9));
t97 = sin(qJ(5));
t147 = t97 * t93;
t109 = t140 - t147;
t45 = t100 * t93 + t97 * t95;
t91 = qJD(1) + qJD(3);
t37 = t45 * t91;
t144 = g(2) * t74 - g(3) * t75;
t101 = cos(qJ(3));
t98 = sin(qJ(3));
t155 = -t159 * t101 - t160 * t98;
t87 = qJDD(1) + qJDD(3);
t12 = t87 * qJ(4) + t91 * qJD(4) - t155;
t79 = t95 * qJDD(2);
t8 = -t93 * t12 + t79;
t9 = t93 * qJDD(2) + t95 * t12;
t157 = -t8 * t93 + t9 * t95;
t143 = t101 * t154 + t98 * t76;
t142 = t93 ^ 2 + t95 ^ 2;
t156 = t142 * t91;
t135 = qJD(1) * t154;
t32 = t101 * t62 - t98 * t135;
t119 = qJD(4) - t32;
t152 = t87 * pkin(3);
t151 = t95 * pkin(4);
t33 = t101 * t135 + t98 * t62;
t150 = t33 * t91;
t38 = t143 * qJD(3);
t149 = t38 * t91;
t148 = t95 * t87;
t145 = t158 * t95;
t139 = t101 * t76;
t137 = t91 * t147;
t133 = t91 * t140;
t136 = qJD(5) * t133 + t45 * t87;
t77 = -pkin(3) - t151;
t131 = -t74 * pkin(3) + t75 * qJ(4);
t130 = t142 * t87;
t126 = t160 * t101 - t159 * t98;
t114 = qJDD(4) - t126;
t10 = t77 * t87 + t114;
t23 = t77 * t91 + t119;
t42 = t45 * qJD(5);
t90 = pkin(9) + qJ(5);
t83 = cos(t90);
t128 = -t10 * t109 + t158 * t83 + t23 * t42;
t127 = -t98 * t154 + t139;
t124 = t109 * t87;
t40 = -pkin(3) - t127;
t121 = t150 + t152;
t102 = cos(qJ(1));
t99 = sin(qJ(1));
t120 = g(2) * t102 + g(3) * t99;
t118 = -t75 * pkin(3) - t74 * qJ(4);
t28 = t91 * qJ(4) + t33;
t21 = t95 * qJD(2) - t93 * t28;
t22 = t93 * qJD(2) + t95 * t28;
t117 = t21 * t93 - t22 * t95;
t116 = -t40 * t87 - t149;
t115 = t144 + t157;
t39 = qJ(4) + t143;
t29 = (-pkin(7) - t39) * t93;
t85 = t95 * pkin(7);
t30 = t95 * t39 + t85;
t113 = t100 * t30 + t97 * t29;
t112 = t100 * t29 - t97 * t30;
t53 = (-pkin(7) - qJ(4)) * t93;
t54 = t95 * qJ(4) + t85;
t111 = t100 * t54 + t97 * t53;
t110 = t100 * t53 - t97 * t54;
t108 = qJD(3) * t139 - t98 * t134;
t107 = t126 + t158;
t14 = t114 - t152;
t41 = t109 * qJD(5);
t82 = sin(t90);
t105 = t10 * t45 - t158 * t82 + t23 * t41;
t104 = -t144 + t155;
t35 = -t133 + t137;
t34 = qJD(4) + t108;
t31 = t40 - t151;
t27 = -t91 * pkin(3) + t119;
t25 = -t42 * qJD(5) + qJDD(5) * t109;
t24 = t41 * qJD(5) + t45 * qJDD(5);
t20 = t91 * t42 - t124;
t19 = -qJD(5) * t137 + t136;
t13 = t14 * t93;
t4 = pkin(7) * t148 + t9;
t3 = t79 + (-pkin(7) * t87 - t12) * t93;
t2 = t19 * t45 + t37 * t41;
t1 = t109 * t19 - t45 * t20 - t41 * t35 - t37 * t42;
t5 = [qJDD(1), t120, -g(2) * t99 + g(3) * t102, (t120 + (t94 ^ 2 + t96 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t87, t127 * t87 + t107 - t149, -t108 * t91 - t143 * t87 + t104, (t116 - t14) * t95 + t145, t13 + (-t116 - t158) * t93, t130 * t39 + t34 * t156 + t115, t14 * t40 + t27 * t38 - g(2) * (-pkin(2) * cos(t92) - t102 * pkin(1) + t118) - g(3) * (-pkin(2) * sin(t92) - t99 * pkin(1) + t131) + t157 * t39 - t117 * t34, t2, t1, t24, t25, 0, t38 * t35 + t31 * t20 + t112 * qJDD(5) + (-qJD(5) * t113 - t34 * t45) * qJD(5) + t128, t38 * t37 + t31 * t19 - t113 * qJDD(5) + (-qJD(5) * t112 - t109 * t34) * qJD(5) + t105; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t8 * t95 + t9 * t93 - g(1), 0, 0, 0, 0, 0, t25, -t24; 0, 0, 0, 0, t87, t107 + t150, t32 * t91 + t104, (t121 - t14) * t95 + t145, t13 + (-t121 - t158) * t93, qJ(4) * t130 + t119 * t156 + t115, -t14 * pkin(3) - t27 * t33 - g(2) * t118 - g(3) * t131 + (t9 * qJ(4) + t119 * t22) * t95 + (-t8 * qJ(4) - t119 * t21) * t93, t2, t1, t24, t25, 0, t77 * t20 + t110 * qJDD(5) - t33 * t35 + (-t111 * qJD(5) - t119 * t45) * qJD(5) + t128, t77 * t19 - t111 * qJDD(5) - t33 * t37 + (-t110 * qJD(5) - t119 * t109) * qJD(5) + t105; 0, 0, 0, 0, 0, 0, 0, -t148, t93 * t87, -t142 * t91 ^ 2, t117 * t91 + t14 - t158, 0, 0, 0, 0, 0, 0.2e1 * t37 * qJD(5) - t124, (-t35 - t137) * qJD(5) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, (t35 - t137) * qJD(5) + t136, t124, qJDD(5), -g(1) * t83 + t100 * t3 - t144 * t82 - t23 * t37 - t97 * t4, g(1) * t82 - t100 * t4 - t144 * t83 + t23 * t35 - t97 * t3;];
tau_reg = t5;
