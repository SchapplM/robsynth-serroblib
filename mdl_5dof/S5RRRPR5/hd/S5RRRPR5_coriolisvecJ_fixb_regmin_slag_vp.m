% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:44
% EndTime: 2021-01-15 23:10:55
% DurationCPUTime: 2.02s
% Computational Cost: add. (3960->281), mult. (10499->390), div. (0->0), fcn. (7720->8), ass. (0->171)
t131 = sin(pkin(9));
t195 = cos(pkin(9));
t136 = cos(qJ(3));
t137 = cos(qJ(2));
t186 = qJD(1) * t137;
t175 = t136 * t186;
t133 = sin(qJ(3));
t134 = sin(qJ(2));
t187 = qJD(1) * t134;
t176 = t133 * t187;
t91 = -t175 + t176;
t93 = -t133 * t186 - t136 * t187;
t65 = t131 * t93 - t195 * t91;
t60 = qJD(5) - t65;
t225 = qJD(5) - t60;
t182 = qJD(1) * qJD(2);
t224 = -0.2e1 * t182;
t126 = pkin(2) * t187;
t150 = -t131 * t91 - t195 * t93;
t217 = t93 * pkin(3);
t33 = pkin(4) * t150 - t65 * pkin(8) - t217;
t124 = t136 * pkin(2) + pkin(3);
t170 = t195 * t133;
t87 = pkin(2) * t170 + t131 * t124;
t81 = pkin(8) + t87;
t223 = (qJD(5) * t81 + t126 + t33) * t60;
t218 = pkin(6) + pkin(7);
t111 = t218 * t137;
t106 = qJD(1) * t111;
t94 = t133 * t106;
t110 = t218 * t134;
t104 = qJD(1) * t110;
t203 = qJD(2) * pkin(2);
t98 = -t104 + t203;
t165 = t136 * t98 - t94;
t88 = t93 * qJ(4);
t51 = t165 + t88;
t135 = cos(qJ(5));
t132 = sin(qJ(5));
t184 = qJD(5) * t132;
t128 = qJD(2) + qJD(3);
t173 = t137 * t182;
t68 = qJD(3) * t175 - t128 * t176 + t136 * t173;
t103 = t133 * t137 + t136 * t134;
t74 = t128 * t103;
t69 = t74 * qJD(1);
t40 = t131 * t68 + t195 * t69;
t222 = -t135 * t40 + t60 * t184;
t177 = qJD(2) * t218;
t161 = qJD(1) * t177;
t99 = t134 * t161;
t221 = (qJD(3) * t98 - t99) * t136;
t220 = qJD(1) * t103;
t41 = -t131 * t69 + t195 * t68;
t55 = t132 * t128 + t135 * t150;
t16 = t55 * qJD(5) + t132 * t41;
t105 = t134 * t177;
t107 = t137 * t177;
t154 = t133 * t110 - t136 * t111;
t143 = t154 * qJD(3) + t133 * t105 - t136 * t107;
t102 = t133 * t134 - t136 * t137;
t73 = t128 * t102;
t141 = t73 * qJ(4) - t103 * qJD(4) + t143;
t185 = qJD(3) * t133;
t192 = t136 * t110;
t149 = -qJD(3) * t192 - t136 * t105 - t133 * t107 - t111 * t185;
t29 = -t74 * qJ(4) - t102 * qJD(4) + t149;
t11 = t131 * t141 + t195 * t29;
t146 = -t103 * qJ(4) - t133 * t111 - t192;
t61 = -t102 * qJ(4) - t154;
t36 = t131 * t146 + t195 * t61;
t193 = t136 * t106;
t156 = -t133 * t98 - t193;
t100 = t137 * t161;
t171 = -t136 * t100 + t133 * t99;
t144 = t156 * qJD(3) + t171;
t142 = -t68 * qJ(4) + t93 * qJD(4) + t144;
t166 = -t133 * t100 - t106 * t185;
t18 = -t69 * qJ(4) - t91 * qJD(4) + t166 + t221;
t6 = t131 * t18 - t195 * t142;
t72 = -t131 * t102 + t195 * t103;
t160 = -t36 * t40 + t6 * t72;
t196 = t91 * qJ(4);
t52 = -t156 - t196;
t201 = t131 * t52;
t45 = t128 * pkin(3) + t51;
t23 = t195 * t45 - t201;
t21 = -t128 * pkin(4) - t23;
t125 = -t137 * pkin(2) - pkin(1);
t109 = t125 * qJD(1);
t75 = t91 * pkin(3) + qJD(4) + t109;
t28 = -pkin(4) * t65 - pkin(8) * t150 + t75;
t71 = t195 * t102 + t131 * t103;
t77 = t102 * pkin(3) + t125;
t34 = t71 * pkin(4) - t72 * pkin(8) + t77;
t43 = -t131 * t74 - t195 * t73;
t7 = t131 * t142 + t195 * t18;
t219 = -(qJD(5) * t34 + t11) * t60 - (qJD(5) * t28 + t7) * t71 + t21 * t43 + t160;
t216 = t21 * t65;
t215 = t21 * t72;
t214 = t34 * t40;
t53 = -t135 * t128 + t132 * t150;
t213 = t53 * t60;
t212 = t55 * t60;
t211 = t55 * t150;
t210 = t60 * t150;
t209 = t150 * t53;
t208 = t93 * t91;
t48 = t195 * t52;
t24 = t131 * t45 + t48;
t155 = t133 * t104 - t193;
t147 = t155 + t196;
t204 = pkin(2) * qJD(3);
t205 = -t136 * t104 - t94;
t56 = t88 + t205;
t207 = -t131 * t56 + t195 * t147 + (t131 * t136 + t170) * t204;
t194 = t131 * t133;
t206 = t131 * t147 + t195 * t56 - (t195 * t136 - t194) * t204;
t202 = t109 * t93;
t200 = t132 * t40;
t198 = t132 * t65;
t183 = qJD(5) * t135;
t15 = t128 * t183 + t135 * t41 - t150 * t184;
t197 = t15 * t132;
t139 = qJD(1) ^ 2;
t191 = t137 * t139;
t138 = qJD(2) ^ 2;
t190 = t138 * t134;
t189 = t138 * t137;
t188 = t134 ^ 2 - t137 ^ 2;
t127 = t134 * t203;
t180 = t72 * t184;
t179 = t60 * t183;
t22 = t128 * pkin(8) + t24;
t157 = t132 * t22 - t135 * t28;
t178 = t150 * t157 + t21 * t184;
t174 = -pkin(2) * t128 - t98;
t58 = t69 * pkin(3) + qJD(2) * t126;
t67 = t74 * pkin(3) + t127;
t164 = t135 * t60;
t163 = pkin(1) * t224;
t9 = t132 * t28 + t135 * t22;
t162 = t6 * t132 + t150 * t9 + t21 * t183;
t159 = t150 * t24 + t23 * t65;
t158 = t40 * t72 + t43 * t60;
t153 = t60 * t198 - t222;
t152 = -t150 * t75 - t6;
t151 = t109 * t91 - t166;
t86 = -pkin(2) * t194 + t195 * t124;
t145 = t206 * t60 - t81 * t40 - t216;
t140 = -t75 * t65 - t7;
t122 = -t195 * pkin(3) - pkin(4);
t121 = t131 * pkin(3) + pkin(8);
t80 = -pkin(4) - t86;
t76 = t126 - t217;
t57 = -t91 ^ 2 + t93 ^ 2;
t47 = (-t220 - t93) * t128;
t46 = t91 * t128 + t68;
t42 = -t131 * t73 + t195 * t74;
t35 = t131 * t61 - t195 * t146;
t26 = t195 * t51 - t201;
t25 = t131 * t51 + t48;
t14 = t42 * pkin(4) - t43 * pkin(8) + t67;
t13 = t40 * pkin(4) - t41 * pkin(8) + t58;
t12 = t135 * t13;
t10 = t131 * t29 - t195 * t141;
t4 = t55 * t164 + t197;
t3 = t60 * t164 + t200 - t211;
t2 = t153 + t209;
t1 = (t15 - t213) * t135 + (-t16 - t212) * t132;
t5 = [0, 0, 0, 0.2e1 * t134 * t173, t188 * t224, t189, -t190, 0, -pkin(6) * t189 + t134 * t163, pkin(6) * t190 + t137 * t163, t68 * t103 + t93 * t73, -t68 * t102 - t103 * t69 + t73 * t91 + t93 * t74, -t73 * t128, -t74 * t128, 0, t125 * t69 + t109 * t74 + t143 * t128 + (qJD(1) * t102 + t91) * t127, t125 * t68 - t109 * t73 - t149 * t128 + (-t93 + t220) * t127, -t10 * t128 + t77 * t40 + t75 * t42 + t58 * t71 - t65 * t67, -t11 * t128 + t150 * t67 + t77 * t41 + t75 * t43 + t58 * t72, t10 * t150 + t11 * t65 - t23 * t43 - t24 * t42 + t35 * t41 - t7 * t71 + t160, -t23 * t10 + t24 * t11 + t6 * t35 + t7 * t36 + t58 * t77 + t75 * t67, -t55 * t180 + (t15 * t72 + t43 * t55) * t135, (-t132 * t55 - t135 * t53) * t43 + (-t197 - t135 * t16 + (t132 * t53 - t135 * t55) * qJD(5)) * t72, t135 * t158 + t15 * t71 - t180 * t60 + t55 * t42, -t132 * t158 - t16 * t71 - t179 * t72 - t53 * t42, t40 * t71 + t60 * t42, t10 * t53 + t12 * t71 + t35 * t16 - t157 * t42 + (t14 * t60 + t214 + (-t22 * t71 - t36 * t60 + t215) * qJD(5)) * t135 + t219 * t132, t10 * t55 + t35 * t15 - t9 * t42 + (-(-qJD(5) * t36 + t14) * t60 - t214 - (-qJD(5) * t22 + t13) * t71 - qJD(5) * t215) * t132 + t219 * t135; 0, 0, 0, -t134 * t191, t188 * t139, 0, 0, 0, t139 * pkin(1) * t134, pkin(1) * t191, -t208, t57, t46, t47, 0, -t91 * t126 + t202 - t155 * t128 + (t133 * t174 - t193) * qJD(3) + t171, t93 * t126 + t205 * t128 + (qJD(3) * t174 + t99) * t136 + t151, -t207 * t128 + t65 * t76 + t152, t206 * t128 - t150 * t76 + t140, t150 * t207 - t206 * t65 - t87 * t40 - t86 * t41 + t159, -t206 * t24 - t207 * t23 - t6 * t86 + t7 * t87 - t75 * t76, t4, t1, t3, t2, -t210, t80 * t16 + t207 * t53 + (-t6 - t223) * t135 + t145 * t132 + t178, t132 * t223 + t145 * t135 + t80 * t15 + t207 * t55 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t57, t46, t47, 0, -t128 * t156 + t144 + t202, t128 * t165 + t151 - t221, t25 * t128 - t217 * t65 + t152, t26 * t128 + t150 * t217 + t140, -t25 * t150 - t26 * t65 + (-t131 * t40 - t195 * t41) * pkin(3) + t159, t23 * t25 - t24 * t26 + (t131 * t7 - t195 * t6 + t75 * t93) * pkin(3), t4, t1, t3, t2, -t210, t122 * t16 - t6 * t135 - (-t132 * t26 + t135 * t33) * t60 - t25 * t53 - t21 * t198 + (-t179 - t200) * t121 + t178, t122 * t15 + (t132 * t33 + t135 * t26) * t60 - t25 * t55 - t135 * t216 + t222 * t121 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 * t150 + t40, t65 * t128 + t41, -t150 ^ 2 - t65 ^ 2, t150 * t23 - t24 * t65 + t58, 0, 0, 0, 0, 0, t153 - t209, -t135 * t60 ^ 2 - t200 - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t53, -t53 ^ 2 + t55 ^ 2, t15 + t213, -t16 + t212, t40, -t132 * t7 - t21 * t55 - t225 * t9 + t12, -t132 * t13 - t135 * t7 + t225 * t157 + t21 * t53;];
tauc_reg = t5;
