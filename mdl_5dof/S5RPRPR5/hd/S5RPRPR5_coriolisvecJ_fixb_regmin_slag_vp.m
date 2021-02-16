% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:55:12
% EndTime: 2021-01-15 11:55:21
% DurationCPUTime: 1.60s
% Computational Cost: add. (1814->215), mult. (5123->329), div. (0->0), fcn. (3720->8), ass. (0->141)
t119 = sin(qJ(5));
t121 = cos(qJ(5));
t158 = qJD(5) * t119;
t116 = sin(pkin(8));
t115 = sin(pkin(9));
t117 = cos(pkin(9));
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t92 = t115 * t122 + t117 * t120;
t129 = qJD(1) * t92;
t72 = t116 * t129;
t177 = t121 * t72;
t154 = qJD(1) * qJD(3);
t146 = t122 * t154;
t139 = t116 * t146;
t165 = qJD(1) * t116;
t150 = t120 * t165;
t140 = t115 * t150;
t70 = qJD(3) * t140 - t117 * t139;
t82 = t92 * t116;
t79 = qJD(3) * t82;
t71 = qJD(1) * t79;
t149 = t122 * t165;
t76 = t117 * t149 - t140;
t11 = -qJD(5) * t177 + t119 * t70 - t121 * t71 - t76 * t158;
t118 = cos(pkin(8));
t157 = t118 * qJD(1);
t193 = qJD(3) - t157;
t100 = -qJD(5) - t193;
t36 = t119 * t76 + t177;
t175 = t36 * t100;
t196 = t11 - t175;
t134 = -t119 * t72 + t121 * t76;
t195 = t134 * t36;
t12 = t134 * qJD(5) - t119 * t71 - t121 * t70;
t176 = t134 * t100;
t194 = -t12 - t176;
t192 = t134 ^ 2 - t36 ^ 2;
t86 = pkin(3) * t150 + qJ(2) * t165 + qJD(4);
t52 = t72 * pkin(4) + t86;
t184 = t72 * pkin(7);
t172 = qJ(2) * t122;
t153 = t118 * t172;
t95 = -pkin(2) * t118 - t116 * pkin(6) - pkin(1);
t85 = t95 * qJD(1) + qJD(2);
t54 = -qJ(4) * t150 + qJD(1) * t153 + t120 * t85;
t178 = t117 * t54;
t168 = t118 * t120;
t171 = qJ(4) * t116;
t128 = -qJ(2) * t168 - t122 * t171;
t81 = t122 * t85;
t53 = t128 * qJD(1) + t81;
t44 = pkin(3) * t193 + t53;
t21 = t115 * t44 + t178;
t10 = t21 - t184;
t9 = t10 * t158;
t191 = t52 * t36 + t9;
t111 = t116 ^ 2;
t185 = 0.2e1 * t111;
t189 = qJD(5) + t100;
t159 = qJD(4) * t116;
t124 = t128 * qJD(3) - t120 * t159;
t155 = qJD(1) * qJD(2);
t147 = t118 * t155;
t160 = qJD(3) * t122;
t179 = t122 * t147 + t85 * t160;
t30 = t124 * qJD(1) + t179;
t164 = qJD(2) * t120;
t148 = t118 * t164;
t127 = -t122 * t159 - t148;
t161 = qJD(3) * t120;
t169 = t116 * t120;
t31 = -t85 * t161 + ((qJ(4) * t169 - t153) * qJD(3) + t127) * qJD(1);
t4 = -t115 * t30 + t117 * t31;
t2 = t71 * pkin(7) + t4;
t5 = t115 * t31 + t117 * t30;
t3 = t70 * pkin(7) + t5;
t152 = -t119 * t3 + t121 * t2;
t188 = -t52 * t134 + t152;
t144 = -t95 + t171;
t187 = t144 * t120 - t153;
t112 = t118 ^ 2;
t186 = -0.2e1 * t111;
t183 = t76 * pkin(7);
t182 = pkin(3) * t115;
t163 = qJD(2) * t122;
t174 = t118 * t163 + t95 * t160;
t42 = t124 + t174;
t43 = t187 * qJD(3) + t127;
t16 = t115 * t43 + t117 * t42;
t48 = t115 * t54;
t24 = t117 * t53 - t48;
t173 = qJ(2) * t120;
t58 = -t144 * t122 + (-pkin(3) - t173) * t118;
t26 = t115 * t58 - t117 * t187;
t181 = t92 * qJD(3) - t118 * t129;
t133 = t115 * t120 - t117 * t122;
t180 = t193 * t133;
t84 = pkin(3) * t139 + t116 * t155;
t123 = qJD(1) ^ 2;
t170 = t111 * t123;
t90 = (pkin(3) * t160 + qJD(2)) * t116;
t93 = pkin(3) * t169 + t116 * qJ(2);
t167 = t111 + t112;
t166 = t120 ^ 2 - t122 ^ 2;
t162 = qJD(3) * t116;
t156 = qJD(3) - t193;
t151 = qJ(2) * t161;
t15 = -t115 * t42 + t117 * t43;
t20 = t117 * t44 - t48;
t23 = -t115 * t53 - t178;
t25 = t115 * t187 + t117 * t58;
t143 = t167 * t123;
t142 = qJD(1) * t156;
t141 = pkin(3) * t149;
t138 = qJD(5) * t92 + t181;
t137 = -qJD(5) * t133 - t180;
t8 = pkin(4) * t193 - t183 + t20;
t136 = -t121 * t10 - t119 * t8;
t135 = t116 * t142;
t83 = t133 * t116;
t45 = -t119 * t83 + t121 * t82;
t46 = -t119 * t82 - t121 * t83;
t132 = 0.2e1 * t167 * t155;
t131 = (-t193 + t157) * t162;
t126 = -t193 ^ 2 - t170;
t108 = t117 * pkin(3) + pkin(4);
t75 = t133 * t162;
t61 = t76 * pkin(4) + t141;
t59 = t82 * pkin(4) + t93;
t55 = -t75 * pkin(4) + t90;
t47 = -t70 * pkin(4) + t84;
t22 = -t82 * pkin(7) + t26;
t19 = -t118 * pkin(4) + t83 * pkin(7) + t25;
t18 = t46 * qJD(5) - t119 * t79 - t121 * t75;
t17 = -t45 * qJD(5) + t119 * t75 - t121 * t79;
t14 = t24 - t183;
t13 = t23 + t184;
t7 = t75 * pkin(7) + t16;
t6 = t79 * pkin(7) + t15;
t1 = [0, 0, 0, 0, t132, qJ(2) * t132, t120 * t146 * t186, t166 * t154 * t185, t120 * t131, t122 * t131, 0, -t193 * t148 + ((-t120 * t95 - t153) * t193 + t85 * t168) * qJD(3) + (t185 + t112) * qJD(1) * (qJ(2) * t160 + t164), -(-t118 * t151 + t174) * t193 + t179 * t118 + (t163 * t185 + (t186 - t112) * t151) * qJD(1), -t4 * t118 + t15 * t193 - t93 * t70 + t90 * t72 - t86 * t75 + t84 * t82, t5 * t118 - t16 * t193 - t93 * t71 + t90 * t76 - t86 * t79 - t84 * t83, -t15 * t76 - t16 * t72 + t20 * t79 + t21 * t75 + t25 * t71 + t26 * t70 + t4 * t83 - t5 * t82, t20 * t15 + t21 * t16 + t4 * t25 + t5 * t26 + t84 * t93 + t86 * t90, t11 * t46 + t134 * t17, -t11 * t45 - t46 * t12 - t134 * t18 - t17 * t36, -t17 * t100 - t11 * t118, t18 * t100 + t12 * t118, 0, -(-t119 * t7 + t121 * t6) * t100 - t152 * t118 + t55 * t36 + t59 * t12 + t47 * t45 + t52 * t18 + (-(-t119 * t19 - t121 * t22) * t100 - t136 * t118) * qJD(5), t59 * t11 - t9 * t118 + t52 * t17 + t55 * t134 + t47 * t46 + ((-qJD(5) * t22 + t6) * t100 + t2 * t118) * t119 + ((qJD(5) * t19 + t7) * t100 + (qJD(5) * t8 + t3) * t118) * t121; 0, 0, 0, 0, -t143, -qJ(2) * t143, 0, 0, 0, 0, 0, t126 * t120, t126 * t122, -t72 * t165 - t181 * t193, -t76 * t165 + t180 * t193, -t133 * t71 + t180 * t72 + t181 * t76 + t92 * t70, -t133 * t4 - t86 * t165 - t180 * t21 - t181 * t20 + t5 * t92, 0, 0, 0, 0, 0, -t36 * t165 + (t119 * t137 + t121 * t138) * t100, -t134 * t165 + (-t119 * t138 + t121 * t137) * t100; 0, 0, 0, 0, 0, 0, t122 * t120 * t170, -t166 * t170, -t120 * t135, -t122 * t135, 0, (-t156 * t85 - t147) * t120 + (-t118 * t142 - t170) * t172, t81 * t193 + (t156 * t157 + t170) * t173 - t179, -t141 * t72 - t193 * t23 - t86 * t76 + t4, -t141 * t76 + t193 * t24 + t86 * t72 - t5, (t21 + t23) * t76 + (-t20 + t24) * t72 + (t115 * t70 + t117 * t71) * pkin(3), -t20 * t23 - t21 * t24 + (t115 * t5 + t117 * t4 - t86 * t149) * pkin(3), t195, t192, t196, t194, 0, (-t119 * t14 + t121 * t13) * t100 - t61 * t36 + (-(-t108 * t119 - t121 * t182) * t100 + t136) * qJD(5) + t188, -t121 * t3 - t119 * t2 - (t119 * t13 + t121 * t14) * t100 - t61 * t134 + ((t108 * t121 - t119 * t182) * t100 - t121 * t8) * qJD(5) + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 * t76 - t70, -(qJD(3) + t193) * t72, -t72 ^ 2 - t76 ^ 2, t20 * t76 + t21 * t72 + t84, 0, 0, 0, 0, 0, t12 - t176, t11 + t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t192, t196, t194, 0, t189 * t136 + t188, (t10 * t100 - t2) * t119 + (-t189 * t8 - t3) * t121 + t191;];
tauc_reg = t1;
