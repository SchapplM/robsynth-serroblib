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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:34:07
% EndTime: 2020-01-03 11:34:11
% DurationCPUTime: 0.80s
% Computational Cost: add. (1050->150), mult. (1737->195), div. (0->0), fcn. (1138->16), ass. (0->107)
t92 = qJ(1) + pkin(8);
t84 = qJ(3) + t92;
t74 = sin(t84);
t75 = cos(t84);
t123 = g(2) * t75 + g(3) * t74;
t101 = cos(qJ(3));
t94 = sin(pkin(8));
t157 = pkin(1) * t94;
t96 = cos(pkin(8));
t76 = t96 * pkin(1) + pkin(2);
t62 = t76 * qJD(1);
t162 = qJD(3) * t62 + qJDD(1) * t157;
t134 = qJD(3) * t157;
t163 = -qJD(1) * t134 + t76 * qJDD(1);
t98 = sin(qJ(3));
t126 = t163 * t101 - t162 * t98;
t115 = qJDD(4) - t126;
t87 = qJDD(1) + qJDD(3);
t153 = t87 * pkin(3);
t14 = t115 - t153;
t108 = -t123 - t14;
t100 = cos(qJ(5));
t95 = cos(pkin(9));
t141 = t100 * t95;
t93 = sin(pkin(9));
t97 = sin(qJ(5));
t148 = t97 * t93;
t110 = t141 - t148;
t45 = t100 * t93 + t97 * t95;
t91 = qJD(1) + qJD(3);
t37 = t45 * t91;
t158 = -t162 * t101 - t163 * t98;
t12 = t87 * qJ(4) + t91 * qJD(4) - t158;
t79 = t95 * qJDD(2);
t8 = -t93 * t12 + t79;
t9 = t93 * qJDD(2) + t95 * t12;
t161 = -t8 * t93 + t9 * t95;
t160 = g(2) * t74 - g(3) * t75;
t144 = t101 * t157 + t98 * t76;
t143 = t93 ^ 2 + t95 ^ 2;
t159 = t143 * t91;
t135 = qJD(1) * t157;
t32 = t101 * t62 - t98 * t135;
t119 = qJD(4) - t32;
t152 = t95 * pkin(4);
t33 = t101 * t135 + t98 * t62;
t151 = t33 * t91;
t38 = t144 * qJD(3);
t150 = t38 * t91;
t149 = t95 * t87;
t146 = t75 * pkin(3) + t74 * qJ(4);
t140 = t101 * t76;
t138 = t91 * t148;
t137 = t108 * t93;
t133 = t91 * t141;
t136 = qJD(5) * t133 + t45 * t87;
t77 = -pkin(3) - t152;
t131 = t143 * t87;
t10 = t77 * t87 + t115;
t23 = t77 * t91 + t119;
t41 = t110 * qJD(5);
t90 = pkin(9) + qJ(5);
t82 = sin(t90);
t129 = t10 * t45 + t123 * t82 + t23 * t41;
t128 = t74 * pkin(3) - t75 * qJ(4);
t127 = -t98 * t157 + t140;
t124 = t110 * t87;
t40 = -pkin(3) - t127;
t121 = -t151 - t153;
t102 = cos(qJ(1));
t99 = sin(qJ(1));
t120 = -g(2) * t102 - g(3) * t99;
t28 = t91 * qJ(4) + t33;
t21 = t95 * qJD(2) - t93 * t28;
t22 = t93 * qJD(2) + t95 * t28;
t118 = t21 * t93 - t22 * t95;
t117 = t40 * t87 + t150;
t116 = -t160 + t161;
t39 = qJ(4) + t144;
t29 = (-pkin(7) - t39) * t93;
t85 = t95 * pkin(7);
t30 = t95 * t39 + t85;
t114 = t100 * t30 + t97 * t29;
t113 = t100 * t29 - t97 * t30;
t53 = (-pkin(7) - qJ(4)) * t93;
t54 = t95 * qJ(4) + t85;
t112 = t100 * t54 + t97 * t53;
t111 = t100 * t53 - t97 * t54;
t109 = qJD(3) * t140 - t98 * t134;
t42 = t45 * qJD(5);
t106 = t123 - t126;
t83 = cos(t90);
t105 = -t10 * t110 - t123 * t83 + t23 * t42;
t104 = t160 + t158;
t35 = -t133 + t138;
t34 = qJD(4) + t109;
t31 = t40 - t152;
t27 = -t91 * pkin(3) + t119;
t25 = -t42 * qJD(5) + qJDD(5) * t110;
t24 = t41 * qJD(5) + t45 * qJDD(5);
t20 = t91 * t42 - t124;
t19 = -qJD(5) * t138 + t136;
t4 = pkin(7) * t149 + t9;
t3 = t79 + (-pkin(7) * t87 - t12) * t93;
t2 = t19 * t45 + t37 * t41;
t1 = t110 * t19 - t45 * t20 - t41 * t35 - t37 * t42;
t5 = [qJDD(1), t120, g(2) * t99 - g(3) * t102, (t120 + (t94 ^ 2 + t96 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t87, t127 * t87 - t106 - t150, -t109 * t91 - t144 * t87 + t104, (t108 - t117) * t95, t117 * t93 - t137, t39 * t131 + t34 * t159 + t116, t14 * t40 + t27 * t38 - g(2) * (pkin(2) * cos(t92) + t102 * pkin(1) + t146) - g(3) * (pkin(2) * sin(t92) + t99 * pkin(1) + t128) + t161 * t39 - t118 * t34, t2, t1, t24, t25, 0, t38 * t35 + t31 * t20 + t113 * qJDD(5) + (-t114 * qJD(5) - t45 * t34) * qJD(5) + t105, t38 * t37 + t31 * t19 - t114 * qJDD(5) + (-t113 * qJD(5) - t110 * t34) * qJD(5) + t129; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t8 * t95 + t9 * t93 - g(1), 0, 0, 0, 0, 0, t25, -t24; 0, 0, 0, 0, t87, -t106 + t151, t32 * t91 + t104, (t108 - t121) * t95, t121 * t93 - t137, qJ(4) * t131 + t119 * t159 + t116, -t14 * pkin(3) - t27 * t33 - g(2) * t146 - g(3) * t128 + (t9 * qJ(4) + t119 * t22) * t95 + (-t8 * qJ(4) - t119 * t21) * t93, t2, t1, t24, t25, 0, t77 * t20 + t111 * qJDD(5) - t33 * t35 + (-t112 * qJD(5) - t119 * t45) * qJD(5) + t105, t77 * t19 - t112 * qJDD(5) - t33 * t37 + (-t111 * qJD(5) - t119 * t110) * qJD(5) + t129; 0, 0, 0, 0, 0, 0, 0, -t149, t93 * t87, -t143 * t91 ^ 2, t118 * t91 + qJDD(4) + t106 - t153, 0, 0, 0, 0, 0, 0.2e1 * t37 * qJD(5) - t124, (-t35 - t138) * qJD(5) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, (t35 - t138) * qJD(5) + t136, t124, qJDD(5), -g(1) * t83 + t100 * t3 + t160 * t82 - t23 * t37 - t97 * t4, g(1) * t82 - t100 * t4 + t160 * t83 + t23 * t35 - t97 * t3;];
tau_reg = t5;
