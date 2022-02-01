% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:11
% DurationCPUTime: 1.03s
% Computational Cost: add. (1466->168), mult. (2383->237), div. (0->0), fcn. (1484->8), ass. (0->132)
t170 = pkin(8) + pkin(9);
t90 = sin(qJ(5));
t91 = sin(qJ(4));
t94 = cos(qJ(5));
t95 = cos(qJ(4));
t57 = t90 * t95 + t94 * t91;
t86 = qJD(4) + qJD(5);
t168 = t57 * t86;
t87 = qJD(1) + qJD(2);
t83 = qJD(3) + t87;
t17 = t83 * t168;
t166 = qJD(5) - t86;
t140 = pkin(1) * qJD(2);
t127 = qJD(1) * t140;
t141 = pkin(1) * qJD(1);
t93 = sin(qJ(2));
t132 = t93 * t141;
t169 = qJD(3) * t132 + t93 * t127;
t147 = t94 * t95;
t153 = t90 * t91;
t56 = -t147 + t153;
t103 = t56 * t86;
t96 = cos(qJ(3));
t150 = t93 * t96;
t92 = sin(qJ(3));
t97 = cos(qJ(2));
t110 = t92 * t97 + t150;
t139 = qJD(3) * t92;
t117 = pkin(2) * t139 - t110 * t141;
t105 = t117 * t83;
t118 = t97 * t127;
t129 = t97 * t141;
t63 = t87 * pkin(2) + t129;
t22 = t92 * t118 + t63 * t139 + t169 * t96;
t167 = -t105 - t22;
t143 = t169 * t92;
t21 = (qJD(3) * t63 + t118) * t96 - t143;
t164 = t83 * pkin(3);
t163 = t95 * pkin(4);
t162 = t96 * pkin(2);
t79 = t97 * pkin(1) + pkin(2);
t49 = pkin(1) * t150 + t92 * t79 + pkin(8);
t161 = -pkin(9) - t49;
t77 = t92 * pkin(2) + pkin(8);
t160 = -pkin(9) - t77;
t137 = qJD(4) * t91;
t81 = pkin(4) * t137;
t15 = t83 * t81 + t22;
t75 = t92 * t132;
t41 = t96 * t63 - t75;
t80 = -pkin(3) - t163;
t26 = t80 * t83 - t41;
t159 = t15 * t56 + t168 * t26;
t158 = -t103 * t26 + t15 * t57;
t138 = qJD(3) * t96;
t28 = t79 * t139 + (t110 * qJD(2) + t93 * t138) * pkin(1);
t157 = t28 * t83;
t42 = t96 * t132 + t92 * t63;
t156 = t42 * t83;
t133 = t83 * t147;
t134 = t83 * t153;
t44 = -t133 + t134;
t46 = t57 * t83;
t155 = t46 * t44;
t98 = qJD(4) ^ 2;
t154 = t77 * t98;
t151 = t92 * t93;
t125 = t170 * t83 + t42;
t24 = t125 * t95;
t149 = t94 * t24;
t146 = t98 * t91;
t136 = qJD(4) * t95;
t35 = -t41 - t164;
t145 = t35 * t136 + t22 * t91;
t144 = t81 + t117;
t142 = t91 ^ 2 - t95 ^ 2;
t135 = pkin(4) * t83 * t91;
t131 = pkin(2) * t138;
t128 = t83 * t136;
t23 = t125 * t91;
t20 = qJD(4) * pkin(4) - t23;
t126 = -pkin(4) * t86 - t20;
t124 = qJD(4) * t170;
t123 = -t35 * t83 - t21;
t122 = qJD(4) * t161;
t121 = qJD(4) * t160;
t48 = pkin(1) * t151 - t96 * t79 - pkin(3);
t116 = -t42 + t81;
t115 = (-qJD(2) + t87) * t141;
t114 = (-qJD(1) - t87) * t140;
t113 = pkin(8) * t98 - t156;
t111 = t49 * t98 + t157;
t109 = qJD(4) * (t41 - t164);
t108 = qJD(4) * t125;
t27 = t79 * t138 + (-t93 * t139 + (t96 * t97 - t151) * qJD(2)) * pkin(1);
t107 = qJD(4) * (t48 * t83 - t27);
t4 = -t91 * t108 + t95 * t21;
t5 = -t95 * t108 - t91 * t21;
t106 = -t26 * t46 - t90 * t4 + t94 * t5;
t16 = qJD(5) * t133 + t94 * t128 - t86 * t134;
t51 = t96 * t129 - t75;
t101 = qJD(4) * ((-pkin(3) - t162) * t83 - t131 + t51);
t100 = t26 * t44 + (t166 * t24 - t5) * t90;
t85 = t95 * pkin(9);
t84 = t98 * t95;
t82 = t83 ^ 2;
t74 = t95 * pkin(8) + t85;
t73 = t170 * t91;
t67 = t80 - t162;
t61 = t95 * t124;
t60 = t91 * t124;
t55 = 0.2e1 * t91 * t128;
t54 = t95 * t77 + t85;
t53 = t160 * t91;
t47 = t48 - t163;
t43 = -0.2e1 * t142 * t83 * qJD(4);
t40 = t95 * t121 - t91 * t131;
t39 = t91 * t121 + t95 * t131;
t38 = t95 * t49 + t85;
t37 = t161 * t91;
t32 = t168 * t86;
t31 = t103 * t86;
t29 = t35 * t137;
t25 = t81 + t28;
t14 = t95 * t122 - t91 * t27;
t13 = t91 * t122 + t95 * t27;
t12 = -t44 ^ 2 + t46 ^ 2;
t7 = t46 * t86 - t17;
t6 = t44 * t86 + t16;
t2 = -t103 * t46 + t16 * t57;
t1 = t103 * t44 - t16 * t56 - t168 * t46 - t57 * t17;
t3 = [0, 0, 0, 0, t93 * t114, t97 * t114, 0, -t22 - t157, -t27 * t83 - t21, t55, t43, t84, -t146, 0, t29 + t91 * t107 + (-t111 - t22) * t95, t95 * t107 + t111 * t91 + t145, t2, t1, -t31, -t32, 0, t25 * t44 + t47 * t17 + (-t90 * t13 + t94 * t14 + (-t37 * t90 - t38 * t94) * qJD(5)) * t86 + t159, t25 * t46 + t47 * t16 - (t94 * t13 + t90 * t14 + (t37 * t94 - t38 * t90) * qJD(5)) * t86 + t158; 0, 0, 0, 0, t93 * t115, t97 * t115, 0, t167, t51 * t83 + (-t118 + (-pkin(2) * t83 - t63) * qJD(3)) * t96 + t143, t55, t43, t84, -t146, 0, t29 + t91 * t101 + (-t154 + t167) * t95, (t154 + t105) * t91 + t95 * t101 + t145, t2, t1, -t31, -t32, 0, t67 * t17 + (-t90 * t39 + t94 * t40 + (-t53 * t90 - t54 * t94) * qJD(5)) * t86 + t51 * t168 + t144 * t44 + t159, t67 * t16 - (t94 * t39 + t90 * t40 + (t53 * t94 - t54 * t90) * qJD(5)) * t86 - t51 * t103 + t144 * t46 + t158; 0, 0, 0, 0, 0, 0, 0, -t22 + t156, t41 * t83 - t21, t55, t43, t84, -t146, 0, t29 + t91 * t109 + (-t113 - t22) * t95, t95 * t109 + t113 * t91 + t145, t2, t1, -t31, -t32, 0, t80 * t17 + (t90 * t60 - t94 * t61 + (t73 * t90 - t74 * t94) * qJD(5)) * t86 + t116 * t44 + t41 * t168 + t159, t80 * t16 - (-t94 * t60 - t90 * t61 + (-t73 * t94 - t74 * t90) * qJD(5)) * t86 + t116 * t46 - t41 * t103 + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91 * t82 * t95, t142 * t82, 0, 0, 0, t123 * t91, t123 * t95, t155, t12, t6, t7, 0, -t44 * t135 - (t90 * t23 - t149) * t86 + (t126 * t90 - t149) * qJD(5) + t106, -t46 * t135 + (t126 * qJD(5) - t23 * t86 - t4) * t94 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t12, t6, t7, 0, t106 + t166 * (-t90 * t20 - t149), (-t166 * t20 - t4) * t94 + t100;];
tauc_reg = t3;
