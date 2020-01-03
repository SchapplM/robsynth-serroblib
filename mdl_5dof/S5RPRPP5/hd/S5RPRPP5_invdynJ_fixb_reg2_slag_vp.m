% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:47
% EndTime: 2019-12-31 18:16:48
% DurationCPUTime: 1.10s
% Computational Cost: add. (827->228), mult. (1456->244), div. (0->0), fcn. (670->4), ass. (0->143)
t136 = qJ(5) * qJDD(1);
t176 = qJD(1) * qJD(5) + t136;
t141 = qJ(5) * qJD(1);
t90 = -pkin(1) - pkin(6);
t45 = t90 * qJD(1) + qJD(2);
t87 = cos(qJ(3));
t20 = (t45 + t141) * t87;
t143 = qJD(4) - t20;
t89 = -pkin(3) - pkin(4);
t10 = qJD(3) * t89 + t143;
t85 = sin(qJ(3));
t131 = t89 * t85;
t111 = -qJ(2) + t131;
t68 = t87 * qJ(4);
t27 = t68 + t111;
t175 = qJD(1) * t27;
t169 = t85 * pkin(3);
t126 = qJ(2) + t169;
t38 = t126 - t68;
t23 = qJD(1) * t38;
t174 = qJDD(1) * qJ(2);
t86 = sin(qJ(1));
t88 = cos(qJ(1));
t156 = g(1) * t88 + g(2) * t86;
t83 = t85 ^ 2;
t84 = t87 ^ 2;
t154 = t83 + t84;
t44 = t90 * qJDD(1) + qJDD(2);
t128 = t154 * t44;
t150 = qJ(5) + t90;
t172 = t89 * qJDD(3);
t39 = t154 * qJDD(1);
t135 = qJD(1) * qJD(3);
t124 = qJ(5) * t135;
t171 = t87 * t124 + t176 * t85;
t164 = t87 * t88;
t166 = t86 * t87;
t170 = g(1) * t166 - g(2) * t164 - g(3) * t85;
t152 = qJ(4) * t85;
t104 = t87 * t89 - t152;
t119 = qJD(3) * pkin(3) - qJD(4);
t165 = t87 * t45;
t22 = -t119 - t165;
t35 = t85 * t45;
t82 = qJD(3) * qJ(4);
t24 = t35 + t82;
t31 = qJD(3) * t165;
t33 = t85 * t44;
t79 = qJDD(3) * qJ(4);
t81 = qJD(3) * qJD(4);
t6 = t31 + t33 + t79 + t81;
t148 = qJD(3) * t85;
t30 = t45 * t148;
t34 = t87 * t44;
t129 = t34 - qJDD(4) - t30;
t147 = qJDD(3) * pkin(3);
t7 = -t129 - t147;
t93 = t6 * t85 - t7 * t87 + (t22 * t85 + t24 * t87) * qJD(3);
t155 = t83 - t84;
t66 = t87 * qJDD(1);
t12 = 0.2e1 * t155 * t135 - 0.2e1 * t85 * t66;
t77 = g(1) * t86;
t76 = g(2) * t88;
t168 = t85 * t86;
t167 = t85 * t88;
t91 = qJD(3) ^ 2;
t163 = t90 * t91;
t162 = pkin(3) * t166 + t86 * t152;
t161 = g(1) * t164 + g(2) * t166;
t144 = t87 * qJD(4);
t160 = -qJ(4) * t66 - qJD(1) * t144;
t130 = 0.2e1 * qJD(1) * qJD(2);
t159 = (t130 + t174) * qJ(2);
t158 = t88 * pkin(1) + t86 * qJ(2);
t157 = t77 - t76;
t92 = qJD(1) ^ 2;
t153 = t91 + t92;
t151 = t92 * qJ(2);
t149 = pkin(1) * qJDD(1);
t64 = t85 * t141;
t14 = t64 + t24;
t146 = t14 * qJD(3);
t145 = t23 * qJD(1);
t13 = qJD(5) + t175;
t142 = qJD(5) + t13;
t140 = qJDD(3) * t87;
t139 = qJDD(3) * t90;
t138 = t85 * qJDD(1);
t137 = t85 * qJDD(3);
t133 = t88 * pkin(6) + t158;
t132 = t90 * t86;
t127 = t87 * t135;
t125 = qJ(2) * t135;
t123 = pkin(3) * t168 + t133;
t122 = -t34 + t170;
t121 = -t13 - t175;
t120 = qJD(3) * t150;
t117 = qJDD(2) - t149;
t116 = -g(2) * t167 + g(3) * t87 - t33;
t115 = t85 * t127;
t114 = -t151 - t157;
t110 = pkin(3) * t87 + t152;
t108 = -qJDD(4) - t122;
t69 = t88 * qJ(2);
t107 = pkin(3) * t167 - t88 * t68 + t69;
t106 = 0.2e1 * t23 * qJD(3);
t105 = t130 + 0.2e1 * t174;
t95 = t104 * qJD(3) - qJD(2);
t11 = t95 + t144;
t94 = t111 * qJDD(1) + qJDD(5) - t160;
t2 = t95 * qJD(1) + t94;
t103 = qJD(1) * t11 + qJDD(1) * t27 + t2;
t47 = t85 * t124;
t100 = -t129 + t47 + t172;
t3 = -t176 * t87 + t100;
t4 = t6 + t171;
t102 = t10 * t148 + t4 * t85 - t157 + (t146 - t3) * t87;
t101 = -t90 * t39 + t157;
t99 = t110 * qJD(3) + qJD(2);
t98 = t105 - t163;
t18 = t99 - t144;
t5 = t99 * qJD(1) + t126 * qJDD(1) + t160;
t97 = -qJD(1) * t18 - qJDD(1) * t38 + t163 - t5;
t96 = -g(1) * t168 - t116 + 0.2e1 * t79 + 0.2e1 * t81;
t50 = t87 * t92 * t85;
t49 = t87 * t139;
t46 = -t84 * t92 - t91;
t43 = -qJDD(3) + t50;
t42 = t155 * t92;
t41 = t91 * t87 + t137;
t40 = -t91 * t85 + t140;
t37 = t150 * t87;
t36 = t150 * t85;
t32 = t110 * qJD(1);
t29 = t84 * qJDD(1) - 0.2e1 * t115;
t28 = t83 * qJDD(1) + 0.2e1 * t115;
t26 = t153 * t87 + t137;
t25 = -t153 * t85 + t140;
t21 = t104 * qJD(1);
t19 = t35 + t64;
t17 = t85 * qJD(5) + t87 * t120;
t16 = -t87 * qJD(5) + t85 * t120;
t1 = [0, 0, 0, 0, 0, qJDD(1), t157, t156, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t149 - t157, t105 - t156, -t117 * pkin(1) - g(1) * (-t86 * pkin(1) + t69) - g(2) * t158 + t159, t29, t12, t40, t28, -t41, 0, 0.2e1 * t87 * t125 + t49 + (-t156 + t98) * t85, (-0.2e1 * t125 - t139) * t85 + t98 * t87 - t161, t101 - t128, -g(1) * (t69 + t132) - g(2) * t133 + t90 * t128 + t159, t29, t40, -t12, 0, t41, t28, t49 + t87 * t106 + (-t156 - t97) * t85, t101 - t93, (t106 + t139) * t85 + t97 * t87 + t161, t5 * t38 + t23 * t18 - g(1) * (t132 + t107) - g(2) * (-t86 * t68 + t123) + t93 * t90, t29, -t12, -t40, t28, -t41, 0, t37 * qJDD(3) + (t121 * t87 - t16) * qJD(3) + (-t103 - t156) * t85, t36 * qJDD(3) + t103 * t87 + (t121 * t85 + t17) * qJD(3) + t161, (t36 * t85 + t37 * t87) * qJDD(1) + (-t16 * t87 + t17 * t85 + (t36 * t87 - t37 * t85) * qJD(3)) * qJD(1) + t102, t4 * t36 + t14 * t17 - t3 * t37 + t10 * t16 + t2 * t27 + t13 * t11 - g(1) * (pkin(4) * t167 + t107) - g(2) * (-t88 * qJ(5) + t123) + (-g(1) * t150 - g(2) * (t85 * pkin(4) - t68)) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t92, t114 + t117, 0, 0, 0, 0, 0, 0, t25, -t26, -t39, t128 + t114, 0, 0, 0, 0, 0, 0, t25, -t39, t26, t93 - t145 - t157, 0, 0, 0, 0, 0, 0, t25, t26, t39, t13 * qJD(1) + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t42, t66, -t50, -t138, qJDD(3), -t87 * t151 - t122, (t151 + t77) * t85 + t116, 0, 0, t50, t66, t42, qJDD(3), t138, -t50, 0.2e1 * t147 + (-t23 * t87 - t32 * t85) * qJD(1) + t108, -t110 * qJDD(1) + ((t24 - t82) * t87 + (t119 + t22) * t85) * qJD(1), (-t23 * t85 + t32 * t87) * qJD(1) + t96, t6 * qJ(4) - t7 * pkin(3) - t23 * t32 - t22 * t35 - g(1) * t162 - g(3) * (t68 - t169) + t110 * t76 + (qJD(4) - t165) * t24, t50, t42, -t66, -t50, -t138, qJDD(3), t87 * t136 + t19 * qJD(3) - t30 - t47 - 0.2e1 * t172 + (t142 * t87 + t21 * t85) * qJD(1) + t108, -t20 * qJD(3) + t31 + (t13 * t85 - t21 * t87) * qJD(1) + t96 + t171, -t104 * qJDD(1) + (-t14 + t19 + t82) * t87 * qJD(1), t4 * qJ(4) + t3 * t89 - t10 * t19 - t13 * t21 - g(1) * (pkin(4) * t166 + t162) - g(3) * (t68 + t131) + t143 * t14 - t104 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t66, t46, -t24 * qJD(3) + t87 * t145 + t170 + t7, 0, 0, 0, 0, 0, 0, t43, t46, -t66, -t146 + (-qJD(1) * t142 - t136) * t87 + t100 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t127 - t138, -0.2e1 * t85 * t135 + t66, -t154 * t92, (t10 * t87 - t14 * t85 + t95) * qJD(1) + t94 + t156;];
tau_reg = t1;
