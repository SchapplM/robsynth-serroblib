% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:17
% EndTime: 2019-12-05 15:09:18
% DurationCPUTime: 0.85s
% Computational Cost: add. (921->184), mult. (2051->223), div. (0->0), fcn. (1549->10), ass. (0->128)
t120 = qJD(1) * qJD(3);
t80 = sin(qJ(3));
t113 = t80 * t120;
t82 = cos(qJ(3));
t126 = qJDD(1) * t82;
t165 = -qJDD(1) * t80 - t82 * t120;
t75 = sin(pkin(8));
t77 = cos(pkin(8));
t115 = -t75 * t126 + t165 * t77;
t94 = -t75 * t113 - t115;
t168 = qJDD(3) * pkin(6) + qJD(2) * qJD(4) + t94;
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t101 = t81 * pkin(4) + t79 * qJ(5);
t96 = pkin(3) + t101;
t167 = t96 * qJD(3);
t36 = t82 * t75 + t80 * t77;
t30 = t36 * qJD(1);
t23 = qJD(3) * pkin(6) + t30;
t147 = t79 * t23;
t16 = t81 * qJD(2) - t147;
t166 = qJD(5) - t16;
t164 = t96 * qJDD(3);
t107 = (-t113 + t126) * t77 + t165 * t75;
t100 = pkin(4) * t79 - qJ(5) * t81;
t28 = t100 * qJD(4) - t79 * qJD(5);
t136 = qJD(3) * t28;
t8 = -t107 + t136 - t164;
t163 = -t8 + t164;
t76 = sin(pkin(7));
t78 = cos(pkin(7));
t102 = g(1) * t78 + g(2) * t76;
t72 = pkin(8) + qJ(3);
t66 = sin(t72);
t162 = t102 * t66;
t13 = -qJD(4) * pkin(4) + t166;
t144 = t81 * t23;
t17 = t79 * qJD(2) + t144;
t14 = qJD(4) * qJ(5) + t17;
t143 = t82 * t77;
t145 = t80 * t75;
t35 = -t143 + t145;
t135 = qJD(3) * t35;
t161 = -qJD(3) * t135 + t36 * qJDD(3);
t131 = qJDD(4) * pkin(4);
t160 = qJDD(5) - t131;
t153 = g(3) * t66;
t67 = cos(t72);
t159 = t102 * t67 + t153;
t158 = pkin(3) * t66;
t157 = pkin(6) * t67;
t83 = qJD(4) ^ 2;
t156 = pkin(6) * t83;
t152 = g(3) * t67;
t151 = t76 * t79;
t150 = t76 * t81;
t149 = t78 * t79;
t148 = t78 * t81;
t146 = t79 * t81;
t142 = t28 - t30;
t141 = t67 * pkin(3) + t66 * pkin(6);
t73 = t79 ^ 2;
t74 = t81 ^ 2;
t140 = -t73 + t74;
t139 = t73 + t74;
t138 = qJD(3) * pkin(3);
t137 = pkin(6) * qJDD(4);
t134 = qJD(3) * t79;
t133 = qJD(4) * t79;
t132 = qJDD(3) * pkin(3);
t130 = t30 * qJD(3);
t123 = t73 * qJDD(3);
t122 = t74 * qJDD(3);
t121 = t81 * qJDD(3);
t118 = qJD(3) * qJD(4);
t117 = qJDD(4) * qJ(5);
t116 = t79 * qJDD(2) + t168 * t81;
t55 = qJD(1) * t145;
t114 = -g(1) * t76 + g(2) * t78;
t111 = t16 + t147;
t110 = qJD(4) * t144 - t81 * qJDD(2) + t168 * t79;
t29 = qJD(1) * t143 - t55;
t22 = -t29 - t138;
t109 = t22 - t138;
t108 = t81 * t130 + t29 * t133 + (g(1) * t148 + g(2) * t150) * t66;
t15 = -t29 - t167;
t106 = t15 - t167;
t104 = t118 * t146;
t103 = t152 + t156;
t99 = t13 * t79 + t14 * t81;
t98 = t16 * t79 - t17 * t81;
t32 = t36 * qJD(3);
t97 = -t32 * qJD(3) - t35 * qJDD(3);
t25 = t67 * t150 - t149;
t27 = t67 * t148 + t151;
t95 = g(1) * t27 + g(2) * t25 - t116;
t93 = t36 * t83 - t97;
t12 = -t107 - t132;
t92 = -t103 - t12 + t132;
t91 = 0.2e1 * t135 * qJD(4) - qJDD(4) * t36;
t90 = -t152 + t162;
t24 = t67 * t151 + t148;
t26 = t67 * t149 - t150;
t89 = g(1) * t26 + g(2) * t24 + t79 * t153 - t110;
t88 = t17 * qJD(4) + t89;
t3 = t117 + (qJD(5) - t147) * qJD(4) + t116;
t4 = t110 + t160;
t87 = t3 * t81 + t4 * t79 + (t13 * t81 - t14 * t79) * qJD(4);
t6 = -t23 * t133 + t116;
t86 = t6 * t81 + t110 * t79 + (-t16 * t81 - t17 * t79) * qJD(4);
t85 = -t139 * t29 * qJD(3) - t159 + (t122 + t123) * pkin(6);
t84 = qJD(3) ^ 2;
t68 = t79 * qJDD(3);
t56 = t84 * t146;
t46 = t78 * t157;
t45 = t76 * t157;
t44 = t140 * t84;
t43 = qJDD(4) * t81 - t83 * t79;
t42 = qJDD(4) * t79 + t83 * t81;
t40 = qJDD(2) + t114;
t37 = t100 * qJD(3);
t34 = -0.2e1 * t104 + t122;
t33 = 0.2e1 * t104 + t123;
t21 = t140 * t118 + t79 * t121;
t5 = t161 * t139;
t2 = t91 * t79 - t93 * t81;
t1 = t93 * t79 + t91 * t81;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t75 ^ 2 + t77 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t97, -t161, 0, -t107 * t35 - t135 * t30 - t29 * t32 + t94 * t36 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t5, t12 * t35 + t135 * t98 + t22 * t32 + t86 * t36 - g(3), 0, 0, 0, 0, 0, 0, t2, t5, -t1, -t135 * t99 + t15 * t32 + t8 * t35 + t87 * t36 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, t43, -t42, 0, -t98 * qJD(4) - t110 * t81 + t6 * t79 + t114, 0, 0, 0, 0, 0, 0, t43, 0, t42, t99 * qJD(4) + t3 * t79 - t4 * t81 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t107 + t90 + t130, (t29 + t55) * qJD(3) + t115 + t159, 0, 0, t33, 0.2e1 * t21, t42, t34, t43, 0, (t109 * qJD(4) - t137) * t79 + t92 * t81 + t108, (-t137 + (t109 + t29) * qJD(4)) * t81 + (-t130 - t92 - t162) * t79, t85 + t86, -t12 * pkin(3) - t22 * t30 - g(1) * (-t78 * t158 + t46) - g(2) * (-t76 * t158 + t45) - g(3) * t141 + t98 * t29 + t86 * pkin(6), t33, t42, -0.2e1 * t21, 0, -t43, t34, (t106 * qJD(4) - t137) * t79 + (-t103 - t136 + t163) * t81 + t108, t85 + t87, (t137 + (-t106 - t29) * qJD(4)) * t81 + (-t142 * qJD(3) - t156 + t163 + t90) * t79, -g(1) * t46 - g(2) * t45 - g(3) * (t101 * t67 + t141) - t99 * t29 + t142 * t15 + t87 * pkin(6) + (-t8 + t162) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t44, t68, t56, t121, qJDD(4), -t22 * t134 + t88, (-qJD(3) * t22 + t153) * t81 + t111 * qJD(4) + t95, 0, 0, -t56, t68, t44, qJDD(4), -t121, t56, 0.2e1 * t131 - qJDD(5) + (-t15 * t79 + t37 * t81) * qJD(3) + t88, -t100 * qJDD(3), -t81 * t153 + 0.2e1 * t117 + (t15 * t81 + t37 * t79) * qJD(3) + (0.2e1 * qJD(5) - t111) * qJD(4) - t95, t3 * qJ(5) - t4 * pkin(4) - t15 * t37 - t13 * t17 - g(1) * (-t26 * pkin(4) + t27 * qJ(5)) - g(2) * (-t24 * pkin(4) + t25 * qJ(5)) + t100 * t153 + t166 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t56, t68, -t73 * t84 - t83, -t14 * qJD(4) + t15 * t134 + t160 - t89;];
tau_reg = t7;
