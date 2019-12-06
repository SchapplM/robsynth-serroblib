% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:29
% EndTime: 2019-12-05 18:56:38
% DurationCPUTime: 2.77s
% Computational Cost: add. (3089->269), mult. (7696->409), div. (0->0), fcn. (6150->8), ass. (0->155)
t108 = qJD(2) + qJD(3);
t113 = sin(qJ(3));
t117 = cos(qJ(2));
t202 = cos(qJ(3));
t149 = qJD(1) * t202;
t114 = sin(qJ(2));
t174 = qJD(1) * t114;
t220 = -t113 * t174 + t117 * t149;
t59 = t220 * t108;
t112 = sin(qJ(4));
t170 = qJD(4) * t112;
t219 = (-t112 * t220 + t170) * pkin(3);
t218 = qJD(5) * t112 + t170;
t204 = qJD(4) + qJD(5);
t111 = sin(qJ(5));
t116 = cos(qJ(4));
t115 = cos(qJ(5));
t176 = t112 * t115;
t88 = t111 * t116 + t176;
t215 = t204 * t88;
t217 = -t220 * t88 + t215;
t178 = t111 * t112;
t207 = t115 * t116;
t86 = t178 - t207;
t64 = t204 * t86;
t190 = -t220 * t86 + t64;
t173 = qJD(1) * t117;
t85 = -t113 * t173 - t114 * t149;
t70 = -t116 * t108 - t112 * t85;
t72 = t108 * t112 - t116 * t85;
t136 = t111 * t70 - t115 * t72;
t37 = t111 * t72 + t115 * t70;
t216 = t136 * t37;
t171 = qJD(3) * t113;
t160 = pkin(1) * t171;
t214 = t160 + t219;
t213 = t136 ^ 2 - t37 ^ 2;
t166 = qJD(5) * t115;
t168 = qJD(5) * t111;
t169 = qJD(4) * t116;
t31 = t108 * t169 + t116 * t59 + t85 * t170;
t32 = t72 * qJD(4) + t112 * t59;
t6 = -t111 * t32 + t115 * t31 - t70 * t166 - t72 * t168;
t80 = qJD(4) - t220;
t78 = qJD(5) + t80;
t212 = t37 * t78 + t6;
t161 = pkin(1) * t173;
t53 = -pkin(2) * t220 + pkin(5) * t85 - t161;
t189 = pkin(1) * qJD(2);
t107 = t113 * t189;
t91 = pkin(5) * t108 + t107;
t42 = t112 * t53 + t116 * t91;
t34 = t42 * t168;
t148 = t202 * qJD(2);
t142 = pkin(1) * t148;
t92 = -t108 * pkin(2) - t142;
t54 = t70 * pkin(3) + t92;
t211 = t54 * t37 + t34;
t122 = t136 * qJD(5) - t111 * t31 - t115 * t32;
t210 = -t136 * t78 + t122;
t134 = qJD(3) * t142;
t164 = qJD(1) * qJD(2);
t146 = t114 * t164;
t89 = t113 * t117 + t202 * t114;
t67 = t108 * t89;
t60 = t67 * qJD(1);
t22 = pkin(1) * t146 + pkin(2) * t60 - pkin(5) * t59;
t41 = -t112 * t91 + t116 * t53;
t10 = qJD(4) * t41 + t112 * t22 + t116 * t134;
t182 = t115 * t42;
t23 = pkin(3) * t80 + t41;
t15 = t111 * t23 + t182;
t121 = -t42 * qJD(4) - t112 * t134 + t116 * t22;
t9 = t60 * pkin(3) + t121;
t8 = t115 * t9;
t209 = -t15 * qJD(5) - t111 * t10 + t54 * t136 + t8;
t140 = qJD(3) * t107;
t208 = t112 * t140 + t92 * t169;
t206 = t218 * t111;
t205 = -qJD(5) * t116 - t169;
t132 = -t113 * t114 + t202 * t117;
t144 = qJD(5) * t23 + t10;
t185 = t112 * t60;
t172 = qJD(2) * t114;
t66 = t108 * t132;
t30 = pkin(1) * t172 + pkin(2) * t67 - pkin(5) * t66;
t62 = -pkin(1) * t117 - pkin(2) * t132 - pkin(5) * t89;
t180 = t116 * t62;
t45 = -pkin(3) * t132 + t180;
t203 = -(qJD(5) * t45 + t112 * t30 + t62 * t169) * t78 + t144 * t132 - t62 * t185;
t201 = t116 * pkin(3);
t200 = t60 * t86;
t199 = t60 * t88;
t198 = t60 * t89;
t197 = t70 * t80;
t196 = t72 * t80;
t195 = t78 * t85;
t194 = t80 * t85;
t193 = t85 * t220;
t192 = t92 * t220;
t188 = pkin(1) * qJD(3);
t187 = t112 * t31;
t183 = t112 * t92;
t181 = t116 * t60;
t143 = t116 * t80;
t175 = t114 ^ 2 - t117 ^ 2;
t158 = t80 * t170;
t157 = t89 * t170;
t156 = t89 * t169;
t155 = t112 * t202;
t154 = t116 * t202;
t153 = t202 * t108;
t147 = t202 * qJD(3);
t141 = t80 * t155;
t139 = -t42 * t85 + t208;
t105 = -t202 * pkin(1) - pkin(2);
t61 = -pkin(2) * t85 - pkin(5) * t220;
t138 = -t107 + t219;
t137 = t66 * t80 + t198;
t133 = -t80 * t169 - t185;
t14 = -t111 * t42 + t115 * t23;
t24 = pkin(3) * t32 + t140;
t131 = t14 * t85 + t217 * t54 + t24 * t86;
t130 = -t15 * t85 - t190 * t54 + t24 * t88;
t129 = t92 * t170 - t183 * t220 + t41 * t85;
t125 = -t45 * t60 - (pkin(3) * t67 + t116 * t30 - t218 * t62) * t78;
t119 = qJD(1) ^ 2;
t118 = qJD(2) ^ 2;
t106 = -pkin(2) - t201;
t104 = pkin(1) * t113 + pkin(5);
t93 = t105 - t201;
t79 = t85 * pkin(3);
t75 = t85 * t161;
t74 = t220 * t161;
t58 = t86 * t89;
t57 = t88 * t89;
t56 = t116 * t61;
t52 = pkin(1) * t174 + t61;
t48 = t112 * t61 + t116 * t142;
t47 = -t220 ^ 2 + t85 ^ 2;
t46 = t60 * t132;
t44 = -t85 * t108 - t60;
t35 = t116 * t52 - t79;
t33 = -t112 * t142 + t56 - t79;
t18 = t80 * t143 + t72 * t85 + t185;
t17 = -t80 ^ 2 * t112 - t70 * t85 + t181;
t13 = t72 * t143 + t187;
t12 = t88 * t66 + (t204 * t207 - t206) * t89;
t11 = -t215 * t89 - t86 * t66;
t5 = -t217 * t78 - t37 * t85 - t200;
t4 = -t136 * t85 - t190 * t78 + t199;
t3 = (t31 - t197) * t116 + (-t32 - t196) * t112;
t2 = t136 * t190 + t6 * t88;
t1 = t122 * t88 + t136 * t217 + t190 * t37 - t6 * t86;
t7 = [0, 0, 0, 0.2e1 * t117 * t146, -0.2e1 * t175 * t164, t118 * t117, -t118 * t114, 0, 0, 0, t59 * t89 - t66 * t85, t132 * t59 + t220 * t66 + t67 * t85 - t198, t66 * t108, -t67 * t108, 0, (-t220 * t172 - t117 * t60 + (-t117 * t67 - t132 * t172) * qJD(1)) * pkin(1), (-t85 * t172 - t117 * t59 + (-t117 * t66 + t89 * t172) * qJD(1)) * pkin(1), -t72 * t157 + (t31 * t89 + t66 * t72) * t116, (-t112 * t72 - t116 * t70) * t66 + (-t187 - t116 * t32 + (t112 * t70 - t116 * t72) * qJD(4)) * t89, t116 * t137 - t132 * t31 - t157 * t80 + t67 * t72, -t112 * t137 + t132 * t32 - t156 * t80 - t67 * t70, t67 * t80 - t46, -t121 * t132 + t30 * t143 - t158 * t62 + t60 * t180 + t66 * t183 + t208 * t89 + t41 * t67, t10 * t132 - t42 * t67 + (-qJD(4) * t62 * t80 + t140 * t89 + t92 * t66) * t116 + (-t92 * qJD(4) * t89 - t30 * t80 - t62 * t60) * t112, -t11 * t136 - t58 * t6, -t11 * t37 + t12 * t136 - t122 * t58 - t57 * t6, t11 * t78 - t132 * t6 - t136 * t67 - t58 * t60, -t12 * t78 - t122 * t132 - t37 * t67 - t57 * t60, t67 * t78 - t46, t54 * t12 + t14 * t67 + t24 * t57 - t8 * t132 + (qJD(5) * t132 * t42 - t125) * t115 + t203 * t111 + (t37 * t156 + (-t122 * t89 + t37 * t66) * t112) * pkin(3), t54 * t11 - t15 * t67 - t24 * t58 - t34 * t132 + (t132 * t9 + t125) * t111 + t203 * t115 + (-t136 * t156 + (-t136 * t66 + t6 * t89) * t112) * pkin(3); 0, 0, 0, -t114 * t119 * t117, t175 * t119, 0, 0, 0, 0, 0, t193, t47, 0, t44, 0, -t75 + (t220 * t174 + (-qJD(2) - t108) * t171) * pkin(1), t74 + (t85 * t174 + (-t148 - t153) * qJD(3)) * pkin(1), t13, t3, t18, t17, t194, -t52 * t143 + t105 * t32 + t133 * t104 + (-t141 + (-qJD(2) * t116 + t70) * t113) * t188 + t129, t72 * t160 + t105 * t31 + (qJD(4) * t104 + t52) * t80 * t112 + (-pkin(1) * t147 * t80 - t104 * t60 - t192) * t116 + t139, t2, t1, t4, t5, t195, -t104 * t199 - t93 * t122 + t131 + ((-t111 * t154 - t115 * t155) * t188 + t104 * t64 - t115 * t35 + t52 * t178) * t78 + t214 * t37, t104 * t200 + t93 * t6 + t130 + (-(-t111 * t155 + t115 * t154) * t188 + t215 * t104 + t111 * t35 + t52 * t176) * t78 - t214 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t47, 0, t44, 0, -t75 + (-qJD(3) + t108) * t107, t74 + (-t147 + t153) * t189, t13, t3, t18, t17, t194, -pkin(2) * t32 - t56 * t80 + t133 * pkin(5) + (t141 + (-qJD(3) * t116 - t70) * t113) * t189 + t129, -t72 * t107 - t116 * t192 - pkin(2) * t31 + t48 * t80 + (t158 - t181) * pkin(5) + t139, t2, t1, t4, t5, t195, -t106 * t122 - (-t111 * t48 + t115 * t33) * t78 + t138 * t37 + ((t205 * t115 + t206) * t78 - t199) * pkin(5) + t131, t106 * t6 + (t111 * t33 + t115 * t48) * t78 - t138 * t136 + (-(t205 * t111 - t112 * t166 - t115 * t170) * t78 + t200) * pkin(5) + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t70, -t70 ^ 2 + t72 ^ 2, t31 + t197, t196 - t32, t60, t42 * t80 - t92 * t72 + t121, t41 * t80 + t92 * t70 - t10, -t216, t213, t212, t210, t60, -(-t111 * t41 - t182) * t78 + (t115 * t60 - t168 * t78 - t72 * t37) * pkin(3) + t209, (-t42 * t78 - t9) * t111 + (t41 * t78 - t144) * t115 + (-t111 * t60 + t136 * t72 - t166 * t78) * pkin(3) + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216, t213, t212, t210, t60, t15 * t78 + t209, -t111 * t9 - t115 * t144 + t14 * t78 + t211;];
tauc_reg = t7;
