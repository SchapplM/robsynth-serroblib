% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:53
% EndTime: 2021-01-15 11:34:00
% DurationCPUTime: 1.42s
% Computational Cost: add. (1594->237), mult. (3137->311), div. (0->0), fcn. (2180->12), ass. (0->137)
t128 = sin(qJ(1));
t131 = cos(qJ(1));
t192 = g(1) * t128 - g(2) * t131;
t203 = qJDD(2) - t192;
t117 = qJD(3) + qJD(5);
t126 = sin(qJ(5));
t129 = cos(qJ(5));
t124 = sin(pkin(8));
t125 = cos(pkin(8));
t127 = sin(qJ(3));
t130 = cos(qJ(3));
t81 = t124 * t130 + t125 * t127;
t71 = t81 * qJD(1);
t179 = t129 * t71;
t168 = qJD(1) * t127;
t157 = t124 * t168;
t167 = qJD(1) * t130;
t75 = t125 * t167 - t157;
t32 = t126 * t75 + t179;
t177 = t32 * t117;
t164 = qJD(5) * t126;
t163 = qJD(1) * qJD(3);
t156 = t127 * t163;
t159 = t130 * qJDD(1);
t155 = t130 * t163;
t160 = t127 * qJDD(1);
t199 = t155 + t160;
t158 = t124 * t159 + t199 * t125;
t39 = t124 * t156 - t158;
t150 = -t124 * t160 + t125 * t159;
t40 = t81 * t163 - t150;
t4 = -qJD(5) * t179 + t126 * t39 - t129 * t40 - t75 * t164;
t202 = t4 + t177;
t145 = -t126 * t71 + t129 * t75;
t201 = t145 * t32;
t178 = t145 * t117;
t5 = t145 * qJD(5) - t126 * t40 - t129 * t39;
t200 = -t5 + t178;
t119 = qJDD(1) * qJ(2);
t120 = qJD(1) * qJD(2);
t152 = g(1) * t131 + g(2) * t128;
t139 = -t152 + 0.2e1 * t120;
t198 = 0.2e1 * t119 + t139;
t175 = qJDD(1) * pkin(1);
t197 = t175 - t203;
t196 = t145 ^ 2 - t32 ^ 2;
t118 = qJ(3) + pkin(8);
t112 = qJ(5) + t118;
t101 = sin(t112);
t102 = cos(t112);
t186 = t71 * pkin(7);
t132 = -pkin(1) - pkin(6);
t91 = t132 * qJD(1) + qJD(2);
t67 = -qJ(4) * t168 + t127 * t91;
t181 = t125 * t67;
t68 = -qJ(4) * t167 + t130 * t91;
t58 = qJD(3) * pkin(3) + t68;
t23 = t124 * t58 + t181;
t12 = t23 - t186;
t84 = pkin(3) * t168 + qJD(1) * qJ(2) + qJD(4);
t45 = t71 * pkin(4) + t84;
t195 = g(3) * t102 + t192 * t101 + t12 * t164 + t45 * t32;
t80 = -t124 * t127 + t125 * t130;
t42 = -t126 * t81 + t129 * t80;
t193 = qJ(4) - t132;
t115 = qJDD(3) + qJDD(5);
t41 = t126 * t80 + t129 * t81;
t165 = qJD(3) * t130;
t166 = qJD(3) * t127;
t73 = t124 * t166 - t125 * t165;
t74 = t81 * qJD(3);
t8 = -t41 * qJD(5) + t126 * t73 - t129 * t74;
t191 = t42 * t115 + t8 * t117;
t190 = qJD(5) - t117;
t162 = qJD(1) * qJD(4);
t90 = t132 * qJDD(1) + qJDD(2);
t82 = t130 * t90;
t21 = -t130 * t162 - t91 * t166 + qJDD(3) * pkin(3) + t82 + (t156 - t159) * qJ(4);
t26 = (-qJ(4) * qJD(1) + t91) * t165 + (-qJ(4) * qJDD(1) - t162 + t90) * t127;
t6 = -t124 * t26 + t125 * t21;
t2 = qJDD(3) * pkin(4) + t40 * pkin(7) + t6;
t7 = t124 * t21 + t125 * t26;
t3 = t39 * pkin(7) + t7;
t189 = g(3) * t101 - t192 * t102 - t126 * t3 + t129 * t2 - t45 * t145;
t185 = t75 * pkin(7);
t184 = pkin(3) * t124;
t183 = g(3) * t127;
t65 = -t130 * qJD(4) + t166 * t193;
t86 = t193 * t130;
t66 = -qJD(3) * t86 - t127 * qJD(4);
t25 = t124 * t65 + t125 * t66;
t54 = t124 * t67;
t30 = t125 * t68 - t54;
t85 = t193 * t127;
t44 = -t124 * t86 - t125 * t85;
t134 = qJD(1) ^ 2;
t174 = t134 * qJ(2);
t173 = t84 * qJD(1);
t123 = t130 ^ 2;
t170 = t127 ^ 2 - t123;
t133 = qJD(3) ^ 2;
t169 = -t133 - t134;
t95 = pkin(3) * t165 + qJD(2);
t161 = qJDD(3) * t127;
t103 = t127 * pkin(3) + qJ(2);
t22 = t125 * t58 - t54;
t24 = -t124 * t66 + t125 * t65;
t29 = -t124 * t68 - t181;
t43 = t124 * t85 - t125 * t86;
t9 = t42 * qJD(5) - t126 * t74 - t129 * t73;
t149 = -t41 * t115 - t9 * t117;
t11 = qJD(3) * pkin(4) - t185 + t22;
t148 = -t126 * t11 - t129 * t12;
t27 = -t80 * pkin(7) + t43;
t28 = -t81 * pkin(7) + t44;
t147 = -t126 * t28 + t129 * t27;
t146 = t126 * t27 + t129 * t28;
t53 = t199 * pkin(3) + qJDD(4) + t119 + t120;
t104 = t125 * pkin(3) + pkin(4);
t144 = t126 * t104 + t129 * t184;
t143 = t129 * t104 - t126 * t184;
t141 = 0.2e1 * qJ(2) * t163 + qJDD(3) * t132;
t140 = -t192 - t174;
t136 = -t22 * t74 - t23 * t73 + t6 * t80 + t7 * t81 - t192;
t135 = -t132 * t133 + t198;
t111 = qJDD(3) * t130;
t110 = cos(t118);
t109 = sin(t118);
t59 = t81 * pkin(4) + t103;
t49 = pkin(3) * t167 + t75 * pkin(4);
t46 = -t73 * pkin(4) + t95;
t17 = -t39 * pkin(4) + t53;
t16 = t30 - t185;
t15 = t29 + t186;
t14 = t73 * pkin(7) + t25;
t13 = t74 * pkin(7) + t24;
t1 = [qJDD(1), t192, t152, -0.2e1 * t175 + t203, t198, t197 * pkin(1) + (t139 + t119) * qJ(2), t123 * qJDD(1) - 0.2e1 * t127 * t155, -0.2e1 * t127 * t159 + 0.2e1 * t170 * t163, -t133 * t127 + t111, -t133 * t130 - t161, 0, t135 * t127 + t141 * t130, -t141 * t127 + t135 * t130, t24 * qJD(3) + t43 * qJDD(3) - t103 * t39 - t152 * t109 + t53 * t81 + t95 * t71 - t84 * t73, -t25 * qJD(3) - t44 * qJDD(3) - t103 * t40 - t110 * t152 + t53 * t80 - t84 * t74 + t95 * t75, -t24 * t75 - t25 * t71 + t44 * t39 + t43 * t40 - t136, t7 * t44 + t23 * t25 + t6 * t43 + t22 * t24 + t53 * t103 + t84 * t95 - g(1) * (t103 * t131 - t128 * t193) - g(2) * (t103 * t128 + t131 * t193), t145 * t8 + t4 * t42, -t145 * t9 - t8 * t32 - t4 * t41 - t42 * t5, t191, t149, 0, t46 * t32 + t59 * t5 + t17 * t41 + t45 * t9 + (-qJD(5) * t146 - t126 * t14 + t129 * t13) * t117 + t147 * t115 - t152 * t101, t46 * t145 + t59 * t4 + t17 * t42 + t45 * t8 - (qJD(5) * t147 + t126 * t13 + t129 * t14) * t117 - t146 * t115 - t152 * t102; 0, 0, 0, qJDD(1), -t134, -t174 - t197, 0, 0, 0, 0, 0, t169 * t127 + t111, t169 * t130 - t161, -qJD(1) * t71 - t74 * qJD(3) + t80 * qJDD(3), -qJD(1) * t75 + t73 * qJD(3) - t81 * qJDD(3), t81 * t39 + t80 * t40 + t73 * t71 + t74 * t75, t136 - t173, 0, 0, 0, 0, 0, -qJD(1) * t32 + t191, -qJD(1) * t145 + t149; 0, 0, 0, 0, 0, 0, t130 * t134 * t127, -t170 * t134, t159, -t160, qJDD(3), t140 * t130 + t183 + t82, g(3) * t130 + (-t140 - t90) * t127, g(3) * t109 - t29 * qJD(3) - t84 * t75 - t192 * t110 + (qJDD(3) * t125 - t71 * t167) * pkin(3) + t6, g(3) * t110 + t30 * qJD(3) + t84 * t71 + t192 * t109 + (-qJDD(3) * t124 - t75 * t167) * pkin(3) - t7, (t23 + t29) * t75 + (-t22 + t30) * t71 + (t124 * t39 + t125 * t40) * pkin(3), -t22 * t29 - t23 * t30 + (t183 + t124 * t7 + t125 * t6 + (-t192 - t173) * t130) * pkin(3), t201, t196, t202, t200, t115, t143 * t115 - t49 * t32 - (-t126 * t16 + t129 * t15) * t117 + (-t117 * t144 + t148) * qJD(5) + t189, -t144 * t115 - t129 * t3 - t126 * t2 - t49 * t145 + (t126 * t15 + t129 * t16) * t117 + (-t129 * t11 - t117 * t143) * qJD(5) + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t75 - t157) * qJD(3) + t158, -0.2e1 * t71 * qJD(3) + t150, -t71 ^ 2 - t75 ^ 2, t22 * t75 + t23 * t71 - t152 + t53, 0, 0, 0, 0, 0, t5 + t178, t4 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t196, t202, t200, t115, t190 * t148 + t189, (-t12 * t117 - t2) * t126 + (-t190 * t11 - t3) * t129 + t195;];
tau_reg = t1;
