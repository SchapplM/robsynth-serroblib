% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR14
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
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR14_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:03
% EndTime: 2021-01-15 12:17:12
% DurationCPUTime: 1.70s
% Computational Cost: add. (1726->266), mult. (3433->353), div. (0->0), fcn. (2320->10), ass. (0->145)
t101 = sin(pkin(8));
t103 = sin(qJ(3));
t106 = cos(qJ(3));
t170 = cos(pkin(8));
t121 = -t101 * t106 - t170 * t103;
t187 = t121 * qJD(1);
t191 = qJD(5) - t187;
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t153 = t105 * qJD(3);
t139 = t170 * t106;
t128 = qJD(1) * t139;
t159 = qJD(1) * t103;
t145 = t101 * t159;
t57 = t128 - t145;
t36 = t102 * t57 - t153;
t198 = t191 * t36;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t190 = g(1) * t104 - g(2) * t107;
t197 = qJDD(2) - t190;
t152 = qJD(1) * qJD(3);
t144 = t103 * t152;
t148 = t106 * qJDD(1);
t151 = qJD(1) * qJD(4);
t157 = qJD(3) * t103;
t108 = -pkin(1) - pkin(6);
t67 = t108 * qJDD(1) + qJDD(2);
t62 = t106 * t67;
t68 = t108 * qJD(1) + qJD(2);
t16 = -t106 * t151 - t68 * t157 + qJDD(3) * pkin(3) + t62 + (t144 - t148) * qJ(4);
t133 = -qJ(4) * qJD(1) + t68;
t156 = qJD(3) * t106;
t22 = t133 * t156 + (-qJ(4) * qJDD(1) - t151 + t67) * t103;
t5 = -t101 * t22 + t170 * t16;
t1 = -qJDD(3) * pkin(4) - t5;
t95 = qJ(3) + pkin(8);
t86 = sin(t95);
t87 = cos(t95);
t113 = g(3) * t86 - t190 * t87;
t158 = qJD(1) * t106;
t79 = t101 * pkin(3) + pkin(7);
t196 = t113 - (pkin(3) * t158 + t57 * pkin(4) - pkin(7) * t187 + qJD(5) * t79) * t191 - t1;
t127 = g(1) * t107 + g(2) * t104;
t97 = qJD(1) * qJD(2);
t117 = -t127 + 0.2e1 * t97;
t96 = qJDD(1) * qJ(2);
t195 = 0.2e1 * t96 + t117;
t135 = t105 * t191;
t132 = qJDD(1) * t170;
t147 = qJD(3) * t128 + t101 * t148 + t103 * t132;
t32 = t101 * t144 - t147;
t30 = -qJDD(5) + t32;
t172 = t102 * t30;
t194 = t135 * t191 - t172;
t169 = qJDD(1) * pkin(1);
t193 = t169 - t197;
t192 = qJD(3) * t187;
t188 = qJ(4) - t108;
t46 = t133 * t103;
t173 = t101 * t46;
t47 = -qJ(4) * t158 + t106 * t68;
t44 = qJD(3) * pkin(3) + t47;
t17 = t170 * t44 - t173;
t13 = -qJD(3) * pkin(4) - t17;
t64 = pkin(3) * t159 + qJD(1) * qJ(2) + qJD(4);
t19 = -pkin(4) * t187 - t57 * pkin(7) + t64;
t6 = t101 * t16 + t170 * t22;
t141 = qJDD(3) * pkin(7) + qJD(5) * t19 + t6;
t116 = -t106 * qJD(4) + t157 * t188;
t134 = t188 * t106;
t45 = -qJD(3) * t134 - t103 * qJD(4);
t21 = t101 * t116 + t170 * t45;
t60 = -t101 * t103 + t139;
t80 = t103 * pkin(3) + qJ(2);
t31 = -pkin(4) * t121 - t60 * pkin(7) + t80;
t65 = t188 * t103;
t35 = -t101 * t134 - t170 * t65;
t136 = qJD(3) * t170;
t56 = -t101 * t156 - t103 * t136;
t185 = t1 * t60 + t13 * t56 - (qJD(5) * t31 + t21) * t191 + t141 * t121 + t35 * t30;
t181 = g(3) * t87;
t180 = g(3) * t103;
t179 = t13 * t60;
t178 = t31 * t30;
t38 = t102 * qJD(3) + t105 * t57;
t177 = t38 * t57;
t176 = t57 * t36;
t155 = qJD(5) * t102;
t149 = t103 * qJDD(1);
t124 = -t101 * t149 + t106 * t132;
t33 = -t124 - t192;
t9 = qJD(5) * t153 + t102 * qJDD(3) - t105 * t33 - t57 * t155;
t175 = t9 * t102;
t42 = t170 * t46;
t18 = t101 * t44 + t42;
t27 = t105 * t30;
t100 = t106 ^ 2;
t171 = -t103 ^ 2 + t100;
t168 = t104 * t102;
t167 = t104 * t105;
t166 = t107 * t102;
t165 = t107 * t105;
t110 = qJD(1) ^ 2;
t164 = t110 * qJ(2);
t163 = t64 * qJD(1);
t109 = qJD(3) ^ 2;
t160 = -t109 - t110;
t72 = pkin(3) * t156 + qJD(2);
t154 = qJD(5) * t105;
t150 = qJDD(3) * t103;
t146 = t60 * t155;
t143 = t106 * t152;
t14 = qJD(3) * pkin(7) + t18;
t41 = qJDD(4) + t96 + t97 + (t143 + t149) * pkin(3);
t8 = -t32 * pkin(4) + t33 * pkin(7) + t41;
t142 = qJD(5) * t14 - t8;
t140 = -t105 * qJDD(3) - t102 * t33;
t131 = -qJD(5) * t121 + qJD(1);
t129 = t56 * t38 + t60 * t9;
t125 = -t191 * t56 + t30 * t60;
t123 = -t27 + (t102 * t187 - t155) * t191;
t122 = -t141 + t181;
t120 = 0.2e1 * qJ(2) * t152 + qJDD(3) * t108;
t118 = -t190 - t164;
t24 = t170 * t47 - t173;
t114 = t79 * t30 + (t13 + t24) * t191;
t55 = t101 * t157 - t106 * t136;
t112 = -t121 * t6 + t17 * t56 - t18 * t55 + t5 * t60 - t190;
t111 = -t108 * t109 + t195;
t90 = qJDD(3) * t106;
t81 = -t170 * pkin(3) - pkin(4);
t52 = t86 * t165 - t168;
t51 = t86 * t166 + t167;
t50 = t86 * t167 + t166;
t49 = -t86 * t168 + t165;
t34 = -t101 * t65 + t170 * t134;
t25 = -t55 * pkin(4) - t56 * pkin(7) + t72;
t23 = t101 * t47 + t42;
t20 = t101 * t45 - t170 * t116;
t10 = t38 * qJD(5) + t140;
t7 = t105 * t8;
t4 = t102 * t19 + t105 * t14;
t3 = -t102 * t14 + t105 * t19;
t2 = [qJDD(1), t190, t127, -0.2e1 * t169 + t197, t195, t193 * pkin(1) + (t117 + t96) * qJ(2), t100 * qJDD(1) - 0.2e1 * t103 * t143, -0.2e1 * t103 * t148 - 0.2e1 * t171 * t152, -t109 * t103 + t90, -t109 * t106 - t150, 0, t111 * t103 + t120 * t106, -t120 * t103 + t111 * t106, -t20 * qJD(3) - t34 * qJDD(3) - t121 * t41 - t127 * t86 - t187 * t72 - t80 * t32 - t64 * t55, -t21 * qJD(3) - t35 * qJDD(3) - t127 * t87 - t80 * t33 + t41 * t60 + t64 * t56 + t72 * t57, t187 * t21 + t20 * t57 + t35 * t32 - t34 * t33 - t112, t6 * t35 + t18 * t21 - t5 * t34 - t17 * t20 + t41 * t80 + t64 * t72 - g(1) * (-t104 * t188 + t80 * t107) - g(2) * (t80 * t104 + t107 * t188), t129 * t105 - t38 * t146, (-t102 * t38 - t105 * t36) * t56 + (-t10 * t105 - t175 + (t102 * t36 - t105 * t38) * qJD(5)) * t60, -t105 * t125 - t121 * t9 - t146 * t191 - t38 * t55, -t154 * t191 * t60 + t10 * t121 + t102 * t125 + t36 * t55, t121 * t30 - t191 * t55, -g(1) * t52 - g(2) * t50 + t34 * t10 + t20 * t36 - t3 * t55 - t7 * t121 + (t25 * t191 - t178 + (t121 * t14 - t191 * t35 + t179) * qJD(5)) * t105 + t185 * t102, g(1) * t51 - g(2) * t49 + t20 * t38 + t34 * t9 + t4 * t55 + (-(-qJD(5) * t35 + t25) * t191 + t178 - t142 * t121 - qJD(5) * t179) * t102 + t185 * t105; 0, 0, 0, qJDD(1), -t110, -t164 - t193, 0, 0, 0, 0, 0, t160 * t103 + t90, t160 * t106 - t150, qJD(1) * t187 + t56 * qJD(3) + t60 * qJDD(3), -qJD(1) * t57 + t55 * qJD(3) + qJDD(3) * t121, -t121 * t32 - t187 * t55 + t60 * t33 - t56 * t57, t112 - t163, 0, 0, 0, 0, 0, -t121 * t172 - t60 * t10 - t56 * t36 + (t102 * t55 - t105 * t131) * t191, -t121 * t27 + (t102 * t131 + t105 * t55) * t191 - t129; 0, 0, 0, 0, 0, 0, t106 * t110 * t103, t171 * t110, t148, -t149, qJDD(3), t118 * t106 + t180 + t62, g(3) * t106 + (-t118 - t67) * t103, t23 * qJD(3) - t64 * t57 + (t170 * qJDD(3) + t158 * t187) * pkin(3) + t113 + t5, t181 + t24 * qJD(3) - t64 * t187 + t190 * t86 + (-qJDD(3) * t101 - t57 * t158) * pkin(3) - t6, (t18 - t23) * t57 - (-t17 + t24) * t187 + (t101 * t32 + t170 * t33) * pkin(3), t17 * t23 - t18 * t24 + (t170 * t5 + t180 + t101 * t6 + (-t190 - t163) * t106) * pkin(3), t38 * t135 + t175, (t9 - t198) * t105 + (-t191 * t38 - t10) * t102, -t177 + t194, t123 + t176, -t191 * t57, t81 * t10 + t114 * t102 + t196 * t105 - t23 * t36 - t3 * t57, -t196 * t102 + t114 * t105 - t23 * t38 + t4 * t57 + t81 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t57 - t145) * qJD(3) + t147, t124 + 0.2e1 * t192, -t187 ^ 2 - t57 ^ 2, t17 * t57 - t18 * t187 - t127 + t41, 0, 0, 0, 0, 0, t123 - t176, -t177 - t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t36, -t36 ^ 2 + t38 ^ 2, t9 + t198, -t140 + (-qJD(5) + t191) * t38, -t30, -g(1) * t49 - g(2) * t51 + t102 * t122 - t13 * t38 - t14 * t154 + t191 * t4 + t7, g(1) * t50 - g(2) * t52 + t102 * t142 + t105 * t122 + t13 * t36 + t191 * t3;];
tau_reg = t2;
