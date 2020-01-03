% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:32
% EndTime: 2019-12-31 20:01:39
% DurationCPUTime: 2.25s
% Computational Cost: add. (3918->344), mult. (10115->443), div. (0->0), fcn. (7056->6), ass. (0->172)
t123 = cos(qJ(2));
t180 = cos(pkin(8));
t157 = t180 * t123;
t109 = qJD(1) * t157;
t119 = sin(pkin(8));
t121 = sin(qJ(2));
t172 = qJD(1) * t121;
t89 = t119 * t172 - t109;
t85 = qJD(4) + t89;
t120 = sin(qJ(4));
t122 = cos(qJ(4));
t170 = qJD(4) * t122;
t171 = qJD(4) * t120;
t197 = -qJ(3) - pkin(6);
t160 = qJD(2) * t197;
t86 = t123 * qJD(3) + t121 * t160;
t80 = t86 * qJD(1);
t129 = -t121 * qJD(3) + t123 * t160;
t81 = t129 * qJD(1);
t42 = t119 * t81 + t180 * t80;
t168 = qJD(1) * qJD(2);
t161 = t121 * t168;
t110 = pkin(2) * t161;
t108 = t119 * t161;
t132 = qJD(2) * t109 - t108;
t101 = t119 * t123 + t180 * t121;
t91 = t101 * qJD(2);
t83 = qJD(1) * t91;
t45 = t83 * pkin(3) - t132 * pkin(7) + t110;
t115 = -t123 * pkin(2) - pkin(1);
t173 = qJD(1) * t115;
t105 = qJD(3) + t173;
t174 = qJD(1) * t101;
t46 = t89 * pkin(3) - pkin(7) * t174 + t105;
t107 = t197 * t123;
t104 = qJD(1) * t107;
t158 = t180 * t104;
t106 = t197 * t121;
t103 = qJD(1) * t106;
t193 = qJD(2) * pkin(2);
t98 = t103 + t193;
t63 = t119 * t98 - t158;
t58 = qJD(2) * pkin(7) + t63;
t156 = t120 * t42 - t122 * t45 + t58 * t170 + t46 * t171;
t21 = t120 * t46 + t122 * t58;
t207 = t21 * t85;
t224 = -t156 + t207;
t73 = t120 * qJD(2) + t122 * t174;
t44 = t73 * qJD(4) + t120 * t132;
t183 = t44 * t122;
t169 = t122 * qJD(2);
t43 = -qJD(4) * t169 - t122 * t132 + t171 * t174;
t184 = t43 * t120;
t71 = t120 * t174 - t169;
t134 = -t119 * t121 + t157;
t94 = t134 * qJD(2);
t223 = t101 * ((t120 * t71 - t122 * t73) * qJD(4) - t183 + t184) - (t120 * t73 + t122 * t71) * t94;
t222 = -0.2e1 * t168;
t20 = -t120 * t58 + t122 * t46;
t4 = t120 * t45 + t122 * t42 + t46 * t170 - t58 * t171;
t221 = -t20 * t85 + t4;
t77 = t120 * t83;
t136 = t85 * t170 + t77;
t78 = t122 * t83;
t219 = -t85 * t171 + t78;
t218 = -t73 * t85 + t44;
t112 = t119 * pkin(2) + pkin(7);
t189 = t112 * t83;
t95 = t119 * t104;
t62 = t180 * t98 + t95;
t57 = -qJD(2) * pkin(3) - t62;
t22 = t71 * pkin(4) - t73 * qJ(5) + t57;
t217 = t22 * t85 - t189;
t187 = t120 * t94;
t216 = t101 * t136 - t134 * t44 + t85 * t187 + t91 * t71;
t61 = -pkin(3) * t134 - t101 * pkin(7) + t115;
t70 = t119 * t106 - t180 * t107;
t194 = t120 * t61 + t122 * t70;
t52 = t119 * t129 + t180 * t86;
t167 = t121 * t193;
t54 = t91 * pkin(3) - t94 * pkin(7) + t167;
t10 = -qJD(4) * t194 - t120 * t52 + t122 * t54;
t214 = t73 ^ 2;
t213 = t85 ^ 2;
t212 = t174 ^ 2;
t211 = t83 * pkin(4);
t210 = pkin(2) * t121;
t16 = t85 * qJ(5) + t21;
t209 = t16 * t85;
t41 = t119 * t80 - t180 * t81;
t69 = -t180 * t106 - t119 * t107;
t206 = t41 * t69;
t205 = t71 * t85;
t204 = t71 * t89;
t203 = t73 * t71;
t201 = t73 * t174;
t200 = t85 * t174;
t199 = t174 * t71;
t198 = t174 * t89;
t147 = pkin(4) * t120 - qJ(5) * t122;
t64 = t119 * t103 - t158;
t196 = t120 * qJD(5) - t85 * t147 + t64;
t195 = -t120 * t44 - t71 * t170;
t166 = pkin(2) * t172;
t53 = pkin(3) * t174 + t89 * pkin(7) + t166;
t65 = t180 * t103 + t95;
t26 = t120 * t53 + t122 * t65;
t192 = t112 * t43;
t191 = t112 * t71;
t190 = t112 * t73;
t188 = t120 * t89;
t186 = t122 * t85;
t185 = t122 * t94;
t182 = t83 * qJ(5);
t181 = t83 * t134;
t125 = qJD(1) ^ 2;
t179 = t123 * t125;
t124 = qJD(2) ^ 2;
t178 = t124 * t121;
t177 = t124 * t123;
t176 = qJD(5) - t20;
t175 = t121 ^ 2 - t123 ^ 2;
t165 = t71 ^ 2 - t214;
t164 = qJD(4) * t112 * t85;
t162 = t121 * t179;
t67 = t73 * t171;
t51 = t119 * t86 - t180 * t129;
t155 = t120 * t85;
t154 = pkin(1) * t222;
t153 = 0.2e1 * t174;
t152 = -t73 * t188 - t67;
t151 = t123 * t161;
t8 = t44 * pkin(4) + t43 * qJ(5) - t73 * qJD(5) + t41;
t150 = -t8 - t164;
t113 = -t180 * pkin(2) - pkin(3);
t149 = t41 + t164;
t148 = t122 * pkin(4) + t120 * qJ(5);
t15 = -t85 * pkin(4) + t176;
t146 = -t120 * t16 + t122 * t15;
t145 = -t120 * t21 - t122 * t20;
t25 = -t120 * t65 + t122 * t53;
t28 = -t120 * t70 + t122 * t61;
t139 = t89 * t186 + t136;
t138 = -t85 * t188 + t219;
t137 = t22 * t73 + t156;
t9 = t120 * t54 + t122 * t52 + t61 * t170 - t70 * t171;
t135 = t85 * t57 - t189;
t133 = (t43 - t204) * t122 + t195;
t131 = t71 * t155 - t183;
t130 = -t138 - t199;
t126 = -t195 * t101 + t71 * t187;
t99 = -t148 + t113;
t88 = t89 ^ 2;
t33 = t73 * pkin(4) + t71 * qJ(5);
t32 = t112 * t183;
t31 = t147 * t101 + t69;
t30 = t85 * t91 - t181;
t24 = pkin(4) * t134 - t28;
t23 = -qJ(5) * t134 + t194;
t19 = -t43 + t205;
t18 = -pkin(4) * t174 - t25;
t17 = qJ(5) * t174 + t26;
t14 = t139 - t201;
t13 = t73 * t186 - t184;
t12 = t147 * t94 + (t148 * qJD(4) - qJD(5) * t122) * t101 + t51;
t11 = t73 * t185 + (-t122 * t43 - t67) * t101;
t7 = -t91 * pkin(4) - t10;
t6 = t91 * qJ(5) - qJD(5) * t134 + t9;
t3 = t101 * t219 + t134 * t43 + t85 * t185 + t73 * t91;
t2 = t156 - t211;
t1 = t85 * qJD(5) + t182 + t4;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t151, t175 * t222, t177, -0.2e1 * t151, -t178, 0, -pkin(6) * t177 + t121 * t154, pkin(6) * t178 + t123 * t154, 0, 0, t132 * t101 + t174 * t94, -t101 * t83 + t132 * t134 - t174 * t91 - t94 * t89, t94 * qJD(2), t89 * t91 - t181, -t91 * qJD(2), 0, t105 * t91 + t115 * t83 + (-t51 + (-qJD(1) * t134 + t89) * t210) * qJD(2), t105 * t94 - t115 * t108 + (t109 * t115 + t153 * t210 - t52) * qJD(2), t41 * t101 + t132 * t69 + t134 * t42 + t174 * t51 - t52 * t89 - t62 * t94 - t63 * t91 - t70 * t83, t206 + t42 * t70 - t62 * t51 + t63 * t52 + (t105 + t173) * t167, t11, t223, t3, t126, -t216, t30, t57 * t187 + t10 * t85 + t156 * t134 + t20 * t91 + t28 * t83 + t69 * t44 + t51 * t71 + (t41 * t120 + t57 * t170) * t101, t57 * t185 + t4 * t134 - t21 * t91 - t194 * t83 - t69 * t43 + t51 * t73 - t9 * t85 + (t41 * t122 - t171 * t57) * t101, -t10 * t73 + t28 * t43 - t194 * t44 - t9 * t71 + t145 * t94 + (-t4 * t120 + t156 * t122 + (t120 * t20 - t122 * t21) * qJD(4)) * t101, t20 * t10 - t156 * t28 + t194 * t4 + t21 * t9 + t57 * t51 + t206, t11, t3, -t223, t30, t216, t126, t22 * t187 + t2 * t134 + t12 * t71 - t15 * t91 - t24 * t83 + t31 * t44 - t7 * t85 + (t8 * t120 + t170 * t22) * t101, -t23 * t44 - t24 * t43 - t6 * t71 + t7 * t73 + t146 * t94 + (-t1 * t120 + t2 * t122 + (-t120 * t15 - t122 * t16) * qJD(4)) * t101, -t22 * t185 - t1 * t134 - t12 * t73 + t16 * t91 + t23 * t83 + t31 * t43 + t6 * t85 + (-t8 * t122 + t171 * t22) * t101, t1 * t23 + t22 * t12 + t15 * t7 + t16 * t6 + t2 * t24 + t8 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t175 * t125, 0, t162, 0, 0, t125 * pkin(1) * t121, pkin(1) * t179, 0, 0, t198, -t88 + t212, -t108 + (t109 + t89) * qJD(2), -t198, 0, 0, t64 * qJD(2) - t105 * t174 - t89 * t166 - t41, t65 * qJD(2) + t105 * t89 - t166 * t174 - t42, (t63 - t64) * t174 + (t65 - t62) * t89 + (-t119 * t83 - t180 * t132) * pkin(2), t62 * t64 - t63 * t65 + (-t105 * t172 + t119 * t42 - t180 * t41) * pkin(2), t13, (-t43 - t204) * t122 + t152 + t195, t14, t131, -t130, -t200, t113 * t44 + t120 * t135 - t122 * t149 - t174 * t20 - t25 * t85 - t64 * t71, -t113 * t43 + t120 * t149 + t122 * t135 + t174 * t21 + t26 * t85 - t64 * t73, t25 * t73 + t26 * t71 - t32 + (-t20 * t89 + t4 + (-t20 + t190) * qJD(4)) * t122 + (-t192 - t21 * t89 + t156 + (-t21 + t191) * qJD(4)) * t120, t41 * t113 - t20 * t25 - t21 * t26 - t57 * t64 + (qJD(4) * t145 + t120 * t156 + t4 * t122) * t112, t13, t14, t67 + (t73 * t89 + t44) * t120 + (t43 + t205) * t122, -t200, t130, t131, t120 * t217 + t150 * t122 + t15 * t174 + t18 * t85 - t196 * t71 + t99 * t44, t17 * t71 - t18 * t73 - t32 + (t15 * t89 + t1 + (t15 + t190) * qJD(4)) * t122 + (-t192 - t16 * t89 + t2 + (-t16 + t191) * qJD(4)) * t120, t150 * t120 - t122 * t217 - t16 * t174 - t17 * t85 + t196 * t73 + t99 * t43, -t15 * t18 - t16 * t17 + t8 * t99 - t196 * t22 + (qJD(4) * t146 + t1 * t122 + t2 * t120) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153 * qJD(2), -t108 + (t109 - t89) * qJD(2), -t88 - t212, t174 * t62 + t63 * t89 + t110, 0, 0, 0, 0, 0, 0, t138 - t199, -t122 * t213 - t201 - t77, t155 * t73 + t133, t221 * t120 + t224 * t122 - t57 * t174, 0, 0, 0, 0, 0, 0, -t120 * t213 - t199 + t78, t133 - t152, t139 + t201, -t22 * t174 + (-t2 + t209) * t122 + (t15 * t85 + t1) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, -t165, t19, -t203, -t218, t83, -t57 * t73 + t224, t57 * t71 - t221, 0, 0, t203, t19, t165, t83, t218, -t203, -t33 * t71 - t137 + t207 + 0.2e1 * t211, pkin(4) * t43 - t44 * qJ(5) + (t16 - t21) * t73 + (t15 - t176) * t71, 0.2e1 * t182 - t22 * t71 + t33 * t73 + (0.2e1 * qJD(5) - t20) * t85 + t4, -t2 * pkin(4) + t1 * qJ(5) - t15 * t21 + t16 * t176 - t22 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t174 + t203, t19, -t213 - t214, t137 - t209 - t211;];
tauc_reg = t5;
