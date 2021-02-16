% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR7
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
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:26
% EndTime: 2021-01-15 12:06:33
% DurationCPUTime: 1.39s
% Computational Cost: add. (1784->252), mult. (3875->344), div. (0->0), fcn. (2712->14), ass. (0->143)
t180 = 2 * qJD(3);
t98 = sin(pkin(8));
t79 = t98 * pkin(1) + pkin(6);
t153 = qJ(4) + t79;
t152 = cos(pkin(9));
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t69 = t79 * qJDD(1);
t113 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + (qJD(3) * qJD(2)) + t69;
t131 = t153 * qJD(1);
t123 = t131 * qJD(3);
t90 = t105 * qJDD(2);
t17 = qJDD(3) * pkin(3) - t102 * t113 - t105 * t123 + t90;
t20 = (qJDD(2) - t123) * t102 + t113 * t105;
t97 = sin(pkin(9));
t3 = t152 * t17 - t97 * t20;
t1 = -qJDD(3) * pkin(4) - t3;
t94 = qJ(1) + pkin(8);
t85 = sin(t94);
t87 = cos(t94);
t129 = g(1) * t87 + g(2) * t85;
t93 = qJ(3) + pkin(9);
t84 = sin(t93);
t86 = cos(t93);
t115 = -g(3) * t86 + t129 * t84;
t150 = qJD(1) * t102;
t136 = t152 * t105;
t75 = qJD(1) * t136;
t57 = t97 * t150 - t75;
t55 = qJD(5) + t57;
t137 = t152 * t102;
t65 = t97 * t105 + t137;
t60 = t65 * qJD(1);
t78 = t97 * pkin(3) + pkin(7);
t179 = t115 - (pkin(3) * t150 + t60 * pkin(4) + t57 * pkin(7) + qJD(5) * t78) * t55 - t1;
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t146 = qJD(1) * qJD(3);
t141 = t102 * t146;
t112 = qJDD(1) * t65 - t97 * t141;
t37 = qJD(3) * t75 + t112;
t45 = t101 * qJD(3) + t104 * t60;
t11 = qJD(5) * t45 - t104 * qJDD(3) + t101 * t37;
t145 = t102 * qJDD(1);
t126 = -qJDD(1) * t136 + t97 * t145;
t59 = t65 * qJD(3);
t36 = qJD(1) * t59 + t126;
t35 = qJDD(5) + t36;
t154 = t97 * t102;
t120 = t136 - t154;
t62 = t120 * qJD(3);
t124 = t35 * t65 + t55 * t62;
t149 = qJD(5) * t101;
t142 = t65 * t149;
t177 = -t104 * t124 + t55 * t142;
t99 = cos(pkin(8));
t173 = t99 * pkin(1);
t83 = t105 * pkin(3) + pkin(2);
t66 = -t83 - t173;
t56 = qJD(1) * t66 + qJD(4);
t25 = t57 * pkin(4) - t60 * pkin(7) + t56;
t4 = t152 * t20 + t97 * t17;
t139 = qJDD(3) * pkin(7) + qJD(5) * t25 + t4;
t48 = t102 * qJD(2) + t105 * t131;
t163 = t97 * t48;
t161 = qJD(3) * pkin(3);
t47 = t105 * qJD(2) - t102 * t131;
t41 = t47 + t161;
t18 = t152 * t41 - t163;
t14 = -qJD(3) * pkin(4) - t18;
t133 = qJD(3) * t153;
t117 = -t102 * qJD(4) - t105 * t133;
t49 = t105 * qJD(4) - t102 * t133;
t24 = t97 * t117 + t152 * t49;
t29 = -pkin(4) * t120 - t65 * pkin(7) + t66;
t63 = t153 * t105;
t33 = t152 * t63 - t153 * t154;
t176 = t1 * t65 - (qJD(5) * t29 + t24) * t55 + t139 * t120 + t14 * t62 - t33 * t35;
t175 = g(3) * t84;
t147 = t104 * qJD(3);
t10 = qJD(5) * t147 + t101 * qJDD(3) + t104 * t37 - t60 * t149;
t172 = -t10 * t120 + t45 * t59;
t171 = pkin(3) * t102;
t170 = g(3) * t105;
t169 = t14 * t65;
t168 = t29 * t35;
t43 = t101 * t60 - t147;
t167 = t43 * t55;
t166 = t45 * t55;
t165 = t45 * t60;
t164 = t60 * t43;
t39 = t152 * t48;
t19 = t97 * t41 + t39;
t95 = t102 ^ 2;
t162 = -t105 ^ 2 + t95;
t160 = t10 * t101;
t159 = t101 * t35;
t158 = t85 * t101;
t157 = t85 * t104;
t156 = t87 * t101;
t155 = t87 * t104;
t81 = -pkin(2) - t173;
t72 = qJD(1) * t81;
t151 = qJDD(2) - g(3);
t148 = qJD(5) * t104;
t144 = t105 * qJDD(1);
t143 = t102 * t161;
t15 = qJD(3) * pkin(7) + t19;
t46 = pkin(3) * t141 + qJDD(1) * t66 + qJDD(4);
t8 = t36 * pkin(4) - t37 * pkin(7) + t46;
t140 = qJD(5) * t15 - t8;
t132 = t104 * t55;
t128 = g(1) * t85 - g(2) * t87;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t127 = g(1) * t103 - g(2) * t106;
t125 = t11 * t120 - t59 * t43;
t122 = t104 * t35 + (-t101 * t57 - t149) * t55;
t121 = -t139 + t175;
t118 = -t72 * qJD(1) + t129 - t69;
t116 = -qJDD(3) * t79 + t72 * t180;
t22 = t152 * t47 - t163;
t114 = -t78 * t35 + (t14 + t22) * t55;
t107 = qJD(3) ^ 2;
t111 = -0.2e1 * qJDD(1) * t81 - t107 * t79 + t128;
t110 = -t65 * t55 * t148 - t101 * t124;
t108 = qJD(1) ^ 2;
t100 = -qJ(4) - pkin(6);
t80 = -t152 * pkin(3) - pkin(4);
t68 = qJDD(3) * t105 - t107 * t102;
t67 = qJDD(3) * t102 + t107 * t105;
t53 = t86 * t155 + t158;
t52 = -t86 * t156 + t157;
t51 = -t86 * t157 + t156;
t50 = t86 * t158 + t155;
t32 = t153 * t137 + t97 * t63;
t27 = t59 * pkin(4) - t62 * pkin(7) + t143;
t23 = -t152 * t117 + t97 * t49;
t21 = t97 * t47 + t39;
t7 = t104 * t8;
t6 = t101 * t25 + t104 * t15;
t5 = -t101 * t15 + t104 * t25;
t2 = [qJDD(1), t127, g(1) * t106 + g(2) * t103, (t127 + (t98 ^ 2 + t99 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t95 * qJDD(1) + 0.2e1 * t105 * t141, 0.2e1 * t102 * t144 - 0.2e1 * t162 * t146, t67, t68, 0, t102 * t116 + t105 * t111, -t102 * t111 + t105 * t116, -t32 * qJDD(3) + t66 * t36 - t46 * t120 + t56 * t59 + t128 * t86 + (t57 * t171 - t23) * qJD(3), -t33 * qJDD(3) + t66 * t37 + t46 * t65 + t56 * t62 - t128 * t84 + (t60 * t171 - t24) * qJD(3), t120 * t4 - t18 * t62 - t19 * t59 + t23 * t60 - t24 * t57 - t3 * t65 + t32 * t37 - t33 * t36 - t129, t4 * t33 + t19 * t24 - t3 * t32 - t18 * t23 + t46 * t66 + t56 * t143 - g(1) * (-t103 * pkin(1) - t87 * t100 - t85 * t83) - g(2) * (t106 * pkin(1) - t85 * t100 + t87 * t83), -t45 * t142 + (t10 * t65 + t45 * t62) * t104, (-t101 * t45 - t104 * t43) * t62 + (-t160 - t104 * t11 + (t101 * t43 - t104 * t45) * qJD(5)) * t65, t172 - t177, t110 + t125, -t120 * t35 + t55 * t59, -g(1) * t51 - g(2) * t53 + t32 * t11 + t23 * t43 + t5 * t59 - t7 * t120 + (t27 * t55 + t168 + (t120 * t15 - t33 * t55 + t169) * qJD(5)) * t104 + t176 * t101, -g(1) * t50 - g(2) * t52 + t32 * t10 + t23 * t45 - t6 * t59 + (-(-qJD(5) * t33 + t27) * t55 - t168 - t140 * t120 - qJD(5) * t169) * t101 + t176 * t104; 0, 0, 0, t151, 0, 0, 0, 0, 0, t68, -t67, -t59 * qJD(3) + qJDD(3) * t120, -t62 * qJD(3) - t65 * qJDD(3), -t120 * t37 - t65 * t36 - t62 * t57 + t59 * t60, t120 * t3 - t18 * t59 + t19 * t62 + t4 * t65 - g(3), 0, 0, 0, 0, 0, t110 - t125, t172 + t177; 0, 0, 0, 0, -t102 * t108 * t105, t162 * t108, t145, t144, qJDD(3), t102 * t118 - t170 + t90, -t151 * t102 + t118 * t105, t21 * qJD(3) - t56 * t60 + (t152 * qJDD(3) - t57 * t150) * pkin(3) + t115 + t3, t175 + t22 * qJD(3) + t56 * t57 + t129 * t86 + (-qJDD(3) * t97 - t60 * t150) * pkin(3) - t4, (t19 - t21) * t60 + (-t18 + t22) * t57 + (-t152 * t37 - t36 * t97) * pkin(3), t18 * t21 - t19 * t22 + (t152 * t3 - t170 + t4 * t97 + (-qJD(1) * t56 + t129) * t102) * pkin(3), t132 * t45 + t160, (t10 - t167) * t104 + (-t11 - t166) * t101, t132 * t55 + t159 - t165, t122 + t164, -t55 * t60, t114 * t101 + t179 * t104 + t80 * t11 - t21 * t43 - t5 * t60, t80 * t10 - t179 * t101 + t114 * t104 - t21 * t45 + t6 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t180 + t126, (t75 - t57) * qJD(3) + t112, -t57 ^ 2 - t60 ^ 2, t18 * t60 + t19 * t57 - t128 + t46, 0, 0, 0, 0, 0, t122 - t164, -t55 ^ 2 * t104 - t159 - t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t43, -t43 ^ 2 + t45 ^ 2, t10 + t167, -t11 + t166, t35, -g(1) * t52 + g(2) * t50 + t101 * t121 - t14 * t45 - t15 * t148 + t6 * t55 + t7, g(1) * t53 - g(2) * t51 + t101 * t140 + t104 * t121 + t14 * t43 + t5 * t55;];
tau_reg = t2;
