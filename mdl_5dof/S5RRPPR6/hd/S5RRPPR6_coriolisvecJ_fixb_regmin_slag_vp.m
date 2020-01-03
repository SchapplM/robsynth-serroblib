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
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 19:33:05
% EndTime: 2019-12-31 19:33:13
% DurationCPUTime: 1.81s
% Computational Cost: add. (2391->256), mult. (6368->376), div. (0->0), fcn. (4696->8), ass. (0->148)
t137 = cos(qJ(2));
t177 = cos(pkin(8));
t155 = t177 * t137;
t122 = qJD(1) * t155;
t132 = sin(pkin(8));
t135 = sin(qJ(2));
t167 = qJD(1) * t135;
t96 = t132 * t167 - t122;
t92 = qJD(5) + t96;
t198 = qJD(5) - t92;
t131 = sin(pkin(9));
t133 = cos(pkin(9));
t112 = t132 * t137 + t177 * t135;
t99 = t112 * qJD(1);
t79 = qJD(2) * t131 + t133 * t99;
t197 = t79 * t96;
t80 = t133 * qJD(2) - t131 * t99;
t196 = t80 * t96;
t136 = cos(qJ(5));
t195 = t136 * t80;
t134 = sin(qJ(5));
t146 = -t134 * t80 - t136 * t79;
t194 = t146 * t92;
t113 = t131 * t136 + t133 * t134;
t178 = t92 * t113;
t164 = qJD(1) * qJD(2);
t193 = -0.2e1 * t164;
t172 = t136 * t133;
t173 = t131 * t134;
t111 = -t172 + t173;
t179 = t92 * t111;
t98 = t112 * qJD(2);
t89 = qJD(1) * t98;
t192 = -t113 * t89 + t179 * t92;
t158 = t135 * t164;
t90 = qJD(2) * t122 - t132 * t158;
t14 = -t146 * qJD(5) + t113 * t90;
t95 = t96 ^ 2;
t191 = pkin(7) * t133;
t36 = t134 * t79 - t195;
t190 = t36 * t99;
t189 = t146 * t99;
t186 = -qJ(3) - pkin(6);
t157 = qJD(2) * t186;
t93 = t137 * qJD(3) + t135 * t157;
t85 = t93 * qJD(1);
t94 = -t135 * qJD(3) + t137 * t157;
t86 = t94 * qJD(1);
t43 = t132 * t85 - t177 * t86;
t119 = t186 * t135;
t120 = t186 * t137;
t75 = -t177 * t119 - t120 * t132;
t188 = t43 * t75;
t124 = pkin(2) * t132 + qJ(4);
t187 = pkin(7) + t124;
t123 = pkin(2) * t158;
t31 = pkin(3) * t89 - qJ(4) * t90 - qJD(4) * t99 + t123;
t44 = t132 * t86 + t177 * t85;
t41 = qJD(2) * qJD(4) + t44;
t9 = t131 * t31 + t133 * t41;
t142 = -t132 * t135 + t155;
t101 = t142 * qJD(2);
t184 = qJD(2) * pkin(2);
t163 = t135 * t184;
t40 = pkin(3) * t98 - qJ(4) * t101 - qJD(4) * t112 + t163;
t57 = t132 * t94 + t177 * t93;
t16 = t131 * t40 + t133 * t57;
t162 = -pkin(2) * t137 - pkin(1);
t151 = t162 * qJD(1);
t118 = qJD(3) + t151;
t47 = t96 * pkin(3) - t99 * qJ(4) + t118;
t115 = qJD(1) * t119;
t109 = t115 + t184;
t116 = qJD(1) * t120;
t156 = t177 * t116;
t67 = t132 * t109 - t156;
t63 = qJD(2) * qJ(4) + t67;
t20 = t131 * t47 + t133 * t63;
t55 = pkin(2) * t167 + pkin(3) * t99 + qJ(4) * t96;
t104 = t132 * t116;
t71 = t177 * t115 + t104;
t25 = t131 * t55 + t133 * t71;
t65 = -pkin(3) * t142 - qJ(4) * t112 + t162;
t76 = t132 * t119 - t177 * t120;
t28 = t131 * t65 + t133 * t76;
t165 = qJD(5) * t136;
t185 = t165 * t80 + t90 * t172;
t182 = t131 * t90;
t181 = t131 * t96;
t180 = t133 * t90;
t176 = t101 * t131;
t175 = t112 * t131;
t174 = t112 * t133;
t139 = qJD(1) ^ 2;
t171 = t137 * t139;
t138 = qJD(2) ^ 2;
t170 = t138 * t135;
t169 = t138 * t137;
t168 = t135 ^ 2 - t137 ^ 2;
t166 = qJD(5) * t112;
t8 = -t131 * t41 + t133 * t31;
t4 = pkin(4) * t89 - pkin(7) * t180 + t8;
t5 = -pkin(7) * t182 + t9;
t161 = -t134 * t5 + t136 * t4;
t19 = -t131 * t63 + t133 * t47;
t160 = -t19 * t96 + t9;
t159 = t20 * t96 + t8;
t15 = -t131 * t57 + t133 * t40;
t24 = -t131 * t71 + t133 * t55;
t27 = -t131 * t76 + t133 * t65;
t56 = t132 * t93 - t177 * t94;
t154 = pkin(1) * t193;
t70 = t115 * t132 - t156;
t153 = -t111 * t89 - t178 * t92;
t127 = -t177 * pkin(2) - pkin(3);
t152 = t134 * t4 + t136 * t5;
t150 = t112 * t43 + t75 * t90;
t12 = pkin(7) * t80 + t20;
t7 = pkin(4) * t96 - pkin(7) * t79 + t19;
t2 = t12 * t136 + t134 * t7;
t149 = t12 * t134 - t136 * t7;
t17 = -pkin(4) * t142 - pkin(7) * t174 + t27;
t21 = -pkin(7) * t175 + t28;
t148 = -t134 * t21 + t136 * t17;
t147 = t134 * t17 + t136 * t21;
t66 = t177 * t109 + t104;
t145 = -qJD(5) * t79 - t182;
t107 = t187 * t131;
t144 = pkin(7) * t181 - qJD(4) * t133 + qJD(5) * t107 + t25;
t108 = t187 * t133;
t143 = pkin(4) * t99 + qJD(4) * t131 + qJD(5) * t108 + t96 * t191 + t24;
t58 = -qJD(2) * pkin(3) + qJD(4) - t66;
t141 = t101 * t58 + t150;
t13 = t145 * t134 + t185;
t140 = -t124 * t89 + t127 * t90 + (-qJD(4) + t58) * t96;
t117 = -t133 * pkin(4) + t127;
t60 = t111 * t112;
t59 = t113 * t112;
t53 = pkin(4) * t175 + t75;
t42 = -pkin(4) * t181 + t70;
t33 = pkin(4) * t176 + t56;
t32 = -pkin(4) * t80 + t58;
t26 = pkin(4) * t182 + t43;
t23 = t113 * t101 + t165 * t174 - t166 * t173;
t22 = -t111 * t101 - t113 * t166;
t10 = -pkin(7) * t176 + t16;
t6 = pkin(4) * t98 - t101 * t191 + t15;
t1 = [0, 0, 0, 0.2e1 * t137 * t158, t168 * t193, t169, -t170, 0, -pkin(6) * t169 + t135 * t154, pkin(6) * t170 + t137 * t154, -t101 * t66 + t142 * t44 + t56 * t99 - t57 * t96 - t67 * t98 - t76 * t89 + t150, t188 + t44 * t76 - t66 * t56 + t67 * t57 + (t118 + t151) * t163, t131 * t141 - t142 * t8 + t15 * t96 + t19 * t98 + t27 * t89 - t56 * t80, t133 * t141 + t142 * t9 - t16 * t96 - t20 * t98 - t28 * t89 + t56 * t79, -t15 * t79 + t16 * t80 + (-t101 * t19 - t112 * t8 - t27 * t90) * t133 + (-t101 * t20 - t112 * t9 - t28 * t90) * t131, t15 * t19 + t16 * t20 + t27 * t8 + t28 * t9 + t56 * t58 + t188, -t13 * t60 - t146 * t22, -t13 * t59 + t14 * t60 + t146 * t23 - t22 * t36, -t13 * t142 - t146 * t98 + t22 * t92 - t60 * t89, t14 * t142 - t23 * t92 - t36 * t98 - t59 * t89, -t142 * t89 + t92 * t98, (-t134 * t10 + t136 * t6) * t92 + t148 * t89 - t161 * t142 - t149 * t98 + t33 * t36 + t53 * t14 + t26 * t59 + t32 * t23 + (t142 * t2 - t147 * t92) * qJD(5), -(t136 * t10 + t134 * t6) * t92 - t147 * t89 + t152 * t142 - t2 * t98 - t33 * t146 + t53 * t13 - t26 * t60 + t32 * t22 + (-t142 * t149 - t148 * t92) * qJD(5); 0, 0, 0, -t135 * t171, t168 * t139, 0, 0, 0, t139 * pkin(1) * t135, pkin(1) * t171, (t67 - t70) * t99 + (-t66 + t71) * t96 + (-t132 * t89 - t177 * t90) * pkin(2), t66 * t70 - t67 * t71 + (-t118 * t167 + t132 * t44 - t177 * t43) * pkin(2), t131 * t140 - t133 * t43 - t19 * t99 - t24 * t96 + t70 * t80, t131 * t43 + t133 * t140 + t20 * t99 + t25 * t96 - t70 * t79, t24 * t79 - t25 * t80 + (qJD(4) * t80 + t160) * t133 + (qJD(4) * t79 - t159) * t131, t127 * t43 - t19 * t24 - t20 * t25 - t58 * t70 + (-t131 * t8 + t133 * t9) * t124 + (-t131 * t19 + t133 * t20) * qJD(4), t13 * t113 + t146 * t179, -t111 * t13 - t113 * t14 + t146 * t178 + t179 * t36, t189 - t192, t153 + t190, -t92 * t99, (-t107 * t136 - t108 * t134) * t89 + t117 * t14 + t26 * t111 + t149 * t99 - t42 * t36 + (t134 * t144 - t136 * t143) * t92 + t178 * t32, -(-t107 * t134 + t108 * t136) * t89 + t117 * t13 + t26 * t113 + t2 * t99 + t42 * t146 + (t134 * t143 + t136 * t144) * t92 - t179 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99 ^ 2 - t95, t66 * t99 + t67 * t96 + t123, -t131 * t95 + t133 * t89 + t80 * t99, -t131 * t89 - t133 * t95 - t79 * t99, (-t180 + t196) * t133 + (-t182 + t197) * t131, t160 * t131 + t159 * t133 - t58 * t99, 0, 0, 0, 0, 0, t153 - t190, t189 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182 + t197, t180 + t196, -t79 ^ 2 - t80 ^ 2, t19 * t79 - t20 * t80 + t43, 0, 0, 0, 0, 0, t14 - t194, t92 * t195 + (-t79 * t92 + t145) * t134 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146 * t36, t146 ^ 2 - t36 ^ 2, t36 * t92 + t13, -t14 - t194, t89, t32 * t146 - t198 * t2 + t161, t198 * t149 + t32 * t36 - t152;];
tauc_reg = t1;
