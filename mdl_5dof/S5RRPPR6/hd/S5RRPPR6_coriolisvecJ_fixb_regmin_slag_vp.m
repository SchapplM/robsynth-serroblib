% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:31
% EndTime: 2021-01-15 19:47:42
% DurationCPUTime: 2.08s
% Computational Cost: add. (2459->271), mult. (6562->395), div. (0->0), fcn. (4816->8), ass. (0->152)
t137 = cos(qJ(2));
t180 = cos(pkin(8));
t156 = t180 * t137;
t121 = qJD(1) * t156;
t132 = sin(pkin(8));
t135 = sin(qJ(2));
t168 = qJD(1) * t135;
t95 = t132 * t168 - t121;
t92 = qJD(5) + t95;
t202 = qJD(5) - t92;
t131 = sin(pkin(9));
t133 = cos(pkin(9));
t111 = t132 * t137 + t135 * t180;
t170 = qJD(1) * t111;
t79 = qJD(2) * t131 + t133 * t170;
t201 = t79 * t95;
t80 = t133 * qJD(2) - t131 * t170;
t200 = t80 * t95;
t136 = cos(qJ(5));
t199 = t136 * t80;
t134 = sin(qJ(5));
t147 = -t134 * t80 - t136 * t79;
t198 = t147 * t92;
t112 = t131 * t136 + t133 * t134;
t181 = t92 * t112;
t165 = qJD(1) * qJD(2);
t197 = -0.2e1 * t165;
t175 = t136 * t133;
t176 = t131 * t134;
t110 = -t175 + t176;
t182 = t92 * t110;
t97 = t111 * qJD(2);
t89 = qJD(1) * t97;
t196 = -t112 * t89 + t182 * t92;
t159 = t135 * t165;
t120 = t132 * t159;
t90 = qJD(2) * t121 - t120;
t14 = -qJD(5) * t147 + t112 * t90;
t94 = t95 ^ 2;
t195 = pkin(2) * t135;
t194 = pkin(7) * t133;
t36 = t134 * t79 - t199;
t193 = t36 * t170;
t192 = t147 * t170;
t189 = -qJ(3) - pkin(6);
t158 = qJD(2) * t189;
t93 = t137 * qJD(3) + t135 * t158;
t85 = t93 * qJD(1);
t141 = -t135 * qJD(3) + t137 * t158;
t86 = t141 * qJD(1);
t43 = t132 * t85 - t180 * t86;
t118 = t189 * t135;
t119 = t189 * t137;
t75 = -t180 * t118 - t119 * t132;
t191 = t43 * t75;
t123 = pkin(2) * t132 + qJ(4);
t190 = pkin(7) + t123;
t122 = pkin(2) * t159;
t31 = pkin(3) * t89 - qJ(4) * t90 - qJD(4) * t170 + t122;
t44 = t132 * t86 + t180 * t85;
t41 = qJD(2) * qJD(4) + t44;
t9 = t131 * t31 + t133 * t41;
t143 = -t132 * t135 + t156;
t100 = t143 * qJD(2);
t187 = qJD(2) * pkin(2);
t163 = t135 * t187;
t40 = pkin(3) * t97 - qJ(4) * t100 - qJD(4) * t111 + t163;
t57 = t132 * t141 + t180 * t93;
t16 = t131 * t40 + t133 * t57;
t127 = -pkin(2) * t137 - pkin(1);
t169 = qJD(1) * t127;
t117 = qJD(3) + t169;
t47 = pkin(3) * t95 - qJ(4) * t170 + t117;
t114 = qJD(1) * t118;
t108 = t114 + t187;
t115 = qJD(1) * t119;
t157 = t180 * t115;
t67 = t132 * t108 - t157;
t63 = qJD(2) * qJ(4) + t67;
t20 = t131 * t47 + t133 * t63;
t164 = pkin(2) * t168;
t55 = pkin(3) * t170 + qJ(4) * t95 + t164;
t103 = t132 * t115;
t71 = t114 * t180 + t103;
t25 = t131 * t55 + t133 * t71;
t65 = -pkin(3) * t143 - qJ(4) * t111 + t127;
t76 = t132 * t118 - t119 * t180;
t28 = t131 * t65 + t133 * t76;
t166 = qJD(5) * t136;
t188 = t166 * t80 + t90 * t175;
t185 = t131 * t90;
t184 = t131 * t95;
t183 = t133 * t90;
t179 = t100 * t131;
t178 = t111 * t131;
t177 = t111 * t133;
t139 = qJD(1) ^ 2;
t174 = t137 * t139;
t138 = qJD(2) ^ 2;
t173 = t138 * t135;
t172 = t138 * t137;
t171 = t135 ^ 2 - t137 ^ 2;
t167 = qJD(5) * t111;
t8 = -t131 * t41 + t133 * t31;
t4 = pkin(4) * t89 - pkin(7) * t183 + t8;
t5 = -pkin(7) * t185 + t9;
t162 = -t134 * t5 + t136 * t4;
t19 = -t131 * t63 + t133 * t47;
t161 = -t19 * t95 + t9;
t160 = t20 * t95 + t8;
t15 = -t131 * t57 + t133 * t40;
t24 = -t131 * t71 + t133 * t55;
t27 = -t131 * t76 + t133 * t65;
t56 = t132 * t93 - t180 * t141;
t155 = pkin(1) * t197;
t154 = 0.2e1 * t170;
t70 = t114 * t132 - t157;
t153 = -t110 * t89 - t181 * t92;
t126 = -pkin(2) * t180 - pkin(3);
t152 = t134 * t4 + t136 * t5;
t151 = t111 * t43 + t75 * t90;
t12 = pkin(7) * t80 + t20;
t7 = pkin(4) * t95 - pkin(7) * t79 + t19;
t2 = t12 * t136 + t134 * t7;
t150 = t12 * t134 - t136 * t7;
t17 = -pkin(4) * t143 - pkin(7) * t177 + t27;
t21 = -pkin(7) * t178 + t28;
t149 = -t134 * t21 + t136 * t17;
t148 = t134 * t17 + t136 * t21;
t66 = t108 * t180 + t103;
t146 = -qJD(5) * t79 - t185;
t106 = t190 * t131;
t145 = pkin(7) * t184 - qJD(4) * t133 + qJD(5) * t106 + t25;
t107 = t190 * t133;
t144 = pkin(4) * t170 + qJD(4) * t131 + qJD(5) * t107 + t95 * t194 + t24;
t58 = -qJD(2) * pkin(3) + qJD(4) - t66;
t142 = t100 * t58 + t151;
t13 = t134 * t146 + t188;
t140 = -t123 * t89 + t126 * t90 + (-qJD(4) + t58) * t95;
t116 = -t133 * pkin(4) + t126;
t60 = t110 * t111;
t59 = t112 * t111;
t53 = pkin(4) * t178 + t75;
t42 = -pkin(4) * t184 + t70;
t33 = pkin(4) * t179 + t56;
t32 = -pkin(4) * t80 + t58;
t26 = pkin(4) * t185 + t43;
t23 = t100 * t112 + t166 * t177 - t167 * t176;
t22 = -t100 * t110 - t112 * t167;
t10 = -pkin(7) * t179 + t16;
t6 = pkin(4) * t97 - t100 * t194 + t15;
t1 = [0, 0, 0, 0.2e1 * t137 * t159, t171 * t197, t172, -t173, 0, -pkin(6) * t172 + t135 * t155, pkin(6) * t173 + t137 * t155, t117 * t97 + t127 * t89 + (-t56 + (-qJD(1) * t143 + t95) * t195) * qJD(2), t100 * t117 + t127 * t90 + (t154 * t195 - t57) * qJD(2), -t100 * t66 + t143 * t44 + t170 * t56 - t57 * t95 - t67 * t97 - t76 * t89 + t151, t191 + t44 * t76 - t56 * t66 + t57 * t67 + (t117 + t169) * t163, t131 * t142 - t143 * t8 + t15 * t95 + t19 * t97 + t27 * t89 - t56 * t80, t133 * t142 + t143 * t9 - t16 * t95 - t20 * t97 - t28 * t89 + t56 * t79, -t15 * t79 + t16 * t80 + (-t100 * t19 - t111 * t8 - t27 * t90) * t133 + (-t100 * t20 - t111 * t9 - t28 * t90) * t131, t15 * t19 + t16 * t20 + t27 * t8 + t28 * t9 + t56 * t58 + t191, -t13 * t60 - t147 * t22, -t13 * t59 + t14 * t60 + t147 * t23 - t22 * t36, -t13 * t143 - t147 * t97 + t22 * t92 - t60 * t89, t14 * t143 - t23 * t92 - t36 * t97 - t59 * t89, -t143 * t89 + t92 * t97, (-t134 * t10 + t136 * t6) * t92 + t149 * t89 - t162 * t143 - t150 * t97 + t33 * t36 + t53 * t14 + t26 * t59 + t32 * t23 + (t143 * t2 - t148 * t92) * qJD(5), -(t136 * t10 + t134 * t6) * t92 - t148 * t89 + t152 * t143 - t2 * t97 - t33 * t147 + t53 * t13 - t26 * t60 + t32 * t22 + (-t143 * t150 - t149 * t92) * qJD(5); 0, 0, 0, -t135 * t174, t171 * t139, 0, 0, 0, t139 * pkin(1) * t135, pkin(1) * t174, qJD(2) * t70 - t117 * t170 - t164 * t95 - t43, qJD(2) * t71 + t117 * t95 - t164 * t170 - t44, (t67 - t70) * t170 + (-t66 + t71) * t95 + (-t132 * t89 - t180 * t90) * pkin(2), t66 * t70 - t67 * t71 + (-t117 * t168 + t132 * t44 - t180 * t43) * pkin(2), t131 * t140 - t133 * t43 - t170 * t19 - t24 * t95 + t70 * t80, t131 * t43 + t133 * t140 + t170 * t20 + t25 * t95 - t70 * t79, t24 * t79 - t25 * t80 + (qJD(4) * t80 + t161) * t133 + (qJD(4) * t79 - t160) * t131, t126 * t43 - t19 * t24 - t20 * t25 - t58 * t70 + (-t131 * t8 + t133 * t9) * t123 + (-t131 * t19 + t133 * t20) * qJD(4), t13 * t112 + t147 * t182, -t110 * t13 - t112 * t14 + t147 * t181 + t182 * t36, t192 - t196, t153 + t193, -t92 * t170, (-t106 * t136 - t107 * t134) * t89 + t116 * t14 + t26 * t110 + t150 * t170 - t42 * t36 + (t134 * t145 - t136 * t144) * t92 + t181 * t32, -(-t106 * t134 + t107 * t136) * t89 + t116 * t13 + t26 * t112 + t2 * t170 + t42 * t147 + (t134 * t144 + t136 * t145) * t92 - t182 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154 * qJD(2), -t120 + (t121 - t95) * qJD(2), -t170 ^ 2 - t94, t170 * t66 + t67 * t95 + t122, -t131 * t94 + t133 * t89 + t170 * t80, -t131 * t89 - t133 * t94 - t170 * t79, (-t183 + t200) * t133 + (-t185 + t201) * t131, t131 * t161 + t133 * t160 - t170 * t58, 0, 0, 0, 0, 0, t153 - t193, t192 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185 + t201, t183 + t200, -t79 ^ 2 - t80 ^ 2, t19 * t79 - t20 * t80 + t43, 0, 0, 0, 0, 0, t14 - t198, t92 * t199 + (-t79 * t92 + t146) * t134 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147 * t36, t147 ^ 2 - t36 ^ 2, t36 * t92 + t13, -t14 - t198, t89, t32 * t147 - t202 * t2 + t162, t202 * t150 + t32 * t36 - t152;];
tauc_reg = t1;
