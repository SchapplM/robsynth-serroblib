% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:06
% EndTime: 2019-12-05 17:26:17
% DurationCPUTime: 3.24s
% Computational Cost: add. (2704->329), mult. (7614->512), div. (0->0), fcn. (6239->12), ass. (0->178)
t111 = sin(pkin(6));
t117 = sin(qJ(3));
t211 = t111 * t117;
t106 = pkin(8) * t211;
t113 = cos(pkin(6));
t121 = cos(qJ(3));
t122 = cos(qJ(2));
t202 = t121 * t122;
t118 = sin(qJ(2));
t206 = t117 * t118;
t134 = -t113 * t206 + t202;
t112 = sin(pkin(5));
t200 = qJD(1) * t112;
t207 = t113 * t121;
t224 = t134 * t200 - (pkin(2) * t207 - t106) * qJD(3);
t151 = pkin(3) * t117 - pkin(9) * t121;
t131 = t151 * qJD(3);
t178 = t118 * t200;
t236 = (-t131 + t178) * t111;
t116 = sin(qJ(4));
t120 = cos(qJ(4));
t196 = qJD(2) * t113;
t164 = qJD(3) + t196;
t198 = qJD(2) * t111;
t177 = t117 * t198;
t235 = -t116 * t177 + t120 * t164;
t67 = qJD(5) - t235;
t195 = qJD(2) * t121;
t176 = t111 * t195;
t104 = -qJD(4) + t176;
t204 = t118 * t121;
t205 = t117 * t122;
t136 = t113 * t204 + t205;
t210 = t111 * t121;
t184 = pkin(8) * t210;
t208 = t113 * t117;
t222 = -t136 * t200 + (pkin(2) * t208 + t184) * qJD(3);
t90 = pkin(8) * t198 + t178;
t221 = qJD(2) * pkin(2);
t97 = t122 * t200 + t221;
t234 = t121 * t90 + t208 * t97;
t188 = qJD(4) * t120;
t190 = qJD(4) * t116;
t77 = t184 + (pkin(2) * t117 + pkin(9)) * t113;
t152 = -pkin(3) * t121 - pkin(9) * t117;
t78 = (-pkin(2) + t152) * t111;
t233 = -t116 * t236 - t120 * t224 + t78 * t188 - t77 * t190;
t223 = t116 * t78 + t120 * t77;
t232 = t223 * qJD(4) - t116 * t224 + t120 * t236;
t231 = t117 * t121;
t114 = cos(pkin(5));
t199 = qJD(1) * t114;
t179 = t111 * t199;
t41 = t117 * t179 + t234;
t32 = pkin(9) * t164 + t41;
t105 = t113 * t199;
t45 = t105 + (qJD(2) * t152 - t97) * t111;
t14 = t116 * t45 + t120 * t32;
t129 = t134 * qJD(2);
t194 = qJD(3) * t111;
t173 = t121 * t194;
t157 = t114 * t173;
t84 = t117 * t90;
t24 = (t207 * t97 - t84) * qJD(3) + (t112 * t129 + t157) * qJD(1);
t54 = (t131 + t178) * t198;
t126 = -qJD(4) * t14 - t116 * t24 + t120 * t54;
t185 = qJD(2) * qJD(3);
t171 = t111 * t185;
t156 = t117 * t171;
t4 = -pkin(4) * t156 - t126;
t72 = t116 * t164 + t120 * t177;
t230 = (pkin(4) * t72 + pkin(10) * t67) * t67 + t4;
t219 = t113 * t97;
t40 = (t179 + t219) * t121 - t84;
t155 = t121 * t171;
t51 = qJD(4) * t72 + t116 * t155;
t115 = sin(qJ(5));
t119 = cos(qJ(5));
t143 = t104 * t115 - t119 * t72;
t50 = qJD(4) * t235 + t120 * t155;
t19 = -qJD(5) * t143 + t115 * t50 - t119 * t156;
t12 = -pkin(10) * t104 + t14;
t31 = -pkin(3) * t164 - t40;
t15 = -pkin(4) * t235 - t72 * pkin(10) + t31;
t148 = t115 * t12 - t119 * t15;
t133 = t116 * t54 + t120 * t24 + t188 * t45 - t190 * t32;
t3 = pkin(10) * t156 + t133;
t130 = t136 * qJD(2);
t193 = qJD(3) * t117;
t174 = t111 * t193;
t158 = t114 * t174;
t25 = t234 * qJD(3) + (t112 * t130 + t158) * qJD(1);
t8 = pkin(4) * t51 - pkin(10) * t50 + t25;
t1 = -qJD(5) * t148 + t115 * t8 + t119 * t3;
t123 = qJD(2) ^ 2;
t47 = t104 * t119 + t115 * t72;
t229 = t47 * t67;
t228 = t143 * t67;
t227 = -pkin(4) * t174 + t232;
t226 = t41 + t104 * (pkin(4) * t116 - pkin(10) * t120);
t79 = t151 * t198;
t225 = t116 * t79 + t120 * t40;
t220 = t104 * t235;
t186 = qJD(5) * t119;
t187 = qJD(5) * t115;
t18 = -t104 * t186 + t115 * t156 + t119 * t50 - t187 * t72;
t218 = t115 * t18;
t217 = t115 * t51;
t216 = t119 * t51;
t215 = t72 * t104;
t214 = t104 * t116;
t213 = t104 * t120;
t108 = t111 ^ 2;
t212 = t108 * t123;
t209 = t112 * t123;
t203 = t120 * t121;
t201 = t117 ^ 2 - t121 ^ 2;
t197 = qJD(2) * t112;
t192 = qJD(3) * t120;
t191 = qJD(4) * t115;
t189 = qJD(4) * t119;
t183 = t67 * t191;
t182 = t67 * t189;
t181 = t115 * t210;
t180 = t118 * t209;
t175 = t118 * t197;
t172 = t111 * t113 * t123;
t167 = t119 * t67;
t103 = -pkin(4) * t120 - pkin(10) * t116 - pkin(3);
t166 = pkin(10) * t177 - qJD(5) * t103 + t225;
t165 = 0.2e1 * t108 * t185;
t163 = qJD(3) + 0.2e1 * t196;
t162 = t108 * t180;
t160 = t111 * t175;
t76 = t106 + (-pkin(2) * t121 - pkin(3)) * t113;
t87 = -t113 * t120 + t116 * t211;
t88 = t113 * t116 + t120 * t211;
t35 = pkin(4) * t87 - pkin(10) * t88 + t76;
t154 = -pkin(10) * t174 - qJD(5) * t35 - t233;
t37 = -pkin(10) * t210 + t223;
t57 = -qJD(4) * t87 + t120 * t173;
t58 = qJD(4) * t88 + t116 * t173;
t153 = -pkin(4) * t58 + pkin(10) * t57 + qJD(5) * t37 - t222;
t6 = t115 * t15 + t119 * t12;
t135 = t113 * t205 + t204;
t56 = t112 * t135 + t114 * t211;
t86 = -t111 * t112 * t122 + t113 * t114;
t34 = t116 * t86 + t120 * t56;
t137 = t113 * t202 - t206;
t55 = -t112 * t137 - t114 * t210;
t147 = t115 * t55 + t119 * t34;
t146 = -t115 * t34 + t119 * t55;
t13 = -t116 * t32 + t120 * t45;
t145 = -t116 * t40 + t120 * t79;
t33 = t116 * t56 - t120 * t86;
t144 = -t116 * t77 + t120 * t78;
t66 = -t111 * t97 + t105;
t141 = -t108 * t221 + t111 * t66;
t140 = -t186 * t67 - t217;
t139 = -t187 * t67 + t216;
t59 = t115 * t88 + t119 * t210;
t11 = pkin(4) * t104 - t13;
t128 = -pkin(10) * t51 + (t11 + t13) * t67;
t127 = qJD(1) * t113 * t175 + qJD(3) * t90;
t2 = -qJD(5) * t6 - t115 * t3 + t119 * t8;
t124 = -t66 * t198 - qJD(3) * t219 + (-t114 * t194 - t122 * t197) * qJD(1);
t65 = (t115 * t117 + t119 * t203) * t198;
t64 = t115 * t120 * t176 - t119 * t177;
t60 = t119 * t88 - t181;
t36 = pkin(4) * t210 - t144;
t30 = t157 + (qJD(3) * t137 + t129) * t112;
t29 = t158 + (qJD(3) * t135 + t130) * t112;
t27 = -qJD(5) * t181 + t115 * t57 - t119 * t174 + t186 * t88;
t26 = -qJD(5) * t59 + t115 * t174 + t119 * t57;
t21 = -pkin(4) * t177 - t145;
t10 = -qJD(4) * t33 + t116 * t160 + t120 * t30;
t9 = qJD(4) * t34 + t116 * t30 - t120 * t160;
t5 = [0, 0, -t180, -t122 * t209, 0, 0, 0, 0, 0, -t121 * t162 + t156 * t86 - t164 * t29, t117 * t162 + t155 * t86 - t164 * t30, 0, 0, 0, 0, 0, t104 * t9 - t156 * t33 - t235 * t29 + t55 * t51, t10 * t104 - t156 * t34 + t29 * t72 + t50 * t55, 0, 0, 0, 0, 0, (-qJD(5) * t147 - t10 * t115 + t119 * t29) * t67 + t146 * t51 + t9 * t47 + t33 * t19, -(qJD(5) * t146 + t10 * t119 + t115 * t29) * t67 - t147 * t51 - t9 * t143 + t33 * t18; 0, 0, 0, 0, t165 * t231, -t201 * t165, t163 * t173, -t163 * t174, 0, (-qJD(2) * t222 - t25) * t113 + (t117 * t141 - t222) * qJD(3), (qJD(2) * t224 - t24) * t113 + (t121 * t141 + t224) * qJD(3), t50 * t88 + t57 * t72, t235 * t57 - t50 * t87 - t51 * t88 - t58 * t72, -t104 * t57 + (-t121 * t50 + (qJD(2) * t88 + t72) * t193) * t111, t104 * t58 + (t121 * t51 + (-qJD(2) * t87 + t235) * t193) * t111, (-t104 * t111 - t108 * t195) * t193, t25 * t87 + t31 * t58 + t76 * t51 - t222 * t235 + t232 * t104 + (-t126 * t121 + (qJD(2) * t144 + t13) * t193) * t111, t25 * t88 + t31 * t57 + t76 * t50 + t222 * t72 + t233 * t104 + (t133 * t121 + (-qJD(2) * t223 - t14) * t193) * t111, -t143 * t26 + t18 * t60, t143 * t27 - t18 * t59 - t19 * t60 - t26 * t47, -t143 * t58 + t18 * t87 + t26 * t67 + t51 * t60, -t19 * t87 - t27 * t67 - t47 * t58 - t51 * t59, t51 * t87 + t58 * t67, (-t115 * t37 + t119 * t35) * t51 + t2 * t87 - t148 * t58 + t36 * t19 + t4 * t59 + t11 * t27 + (t115 * t154 - t119 * t153) * t67 + t227 * t47, -(t115 * t35 + t119 * t37) * t51 - t1 * t87 - t6 * t58 + t36 * t18 + t4 * t60 + t11 * t26 + (t115 * t153 + t119 * t154) * t67 - t227 * t143; 0, 0, 0, 0, -t212 * t231, t201 * t212, -t121 * t172, t117 * t172, 0, t117 * t124 - t121 * t127 + t164 * t41, t117 * t127 + t121 * t124 + t164 * t40, t116 * t50 - t213 * t72, (t50 - t220) * t120 + (-t51 + t215) * t116, -t104 * t188 + (t104 * t203 + (qJD(3) * t116 - t72) * t117) * t198, t104 * t190 + (-t121 * t214 + (-t235 + t192) * t117) * t198, t104 * t177, -pkin(3) * t51 - t25 * t120 + t145 * t104 + t41 * t235 + (pkin(9) * t213 + t116 * t31) * qJD(4) + (-t117 * t13 + (-pkin(9) * t193 - t121 * t31) * t116) * t198, -pkin(3) * t50 + t25 * t116 - t225 * t104 - t41 * t72 + (-pkin(9) * t214 + t120 * t31) * qJD(4) + (-t31 * t203 + (-pkin(9) * t192 + t14) * t117) * t198, t116 * t119 * t18 - (-t116 * t187 + t119 * t188 - t65) * t143, t47 * t65 - t143 * t64 + (t115 * t143 - t119 * t47) * t188 + (-t218 - t119 * t19 + (t115 * t47 + t119 * t143) * qJD(5)) * t116, -t65 * t67 + (-t18 + t182) * t120 + (t104 * t143 + t139) * t116, t64 * t67 + (t19 - t183) * t120 + (t104 * t47 + t140) * t116, -t120 * t51 - t214 * t67, t103 * t216 - t11 * t64 - t21 * t47 + (t115 * t166 - t119 * t226) * t67 + (t11 * t191 - t2 + (qJD(4) * t47 + t140) * pkin(9)) * t120 + (t11 * t186 + t4 * t115 + t104 * t148 + (t19 + t183) * pkin(9)) * t116, -t103 * t217 - t11 * t65 + t21 * t143 + (t115 * t226 + t119 * t166) * t67 + (t11 * t189 + t1 + (-qJD(4) * t143 - t139) * pkin(9)) * t120 + (-t11 * t187 + t4 * t119 + t104 * t6 + (t18 + t182) * pkin(9)) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t235, -t235 ^ 2 + t72 ^ 2, t50 + t220, -t215 - t51, t156, -t104 * t14 - t31 * t72 + t126, -t104 * t13 - t235 * t31 - t133, -t143 * t167 + t218, (t18 - t229) * t119 + (-t19 + t228) * t115, t143 * t72 + t167 * t67 + t217, -t115 * t67 ^ 2 + t47 * t72 + t216, -t67 * t72, -pkin(4) * t19 + t115 * t128 - t119 * t230 - t14 * t47 + t148 * t72, -pkin(4) * t18 + t115 * t230 + t119 * t128 + t14 * t143 + t6 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 * t47, t143 ^ 2 - t47 ^ 2, t18 + t229, -t19 - t228, t51, t11 * t143 + t6 * t67 + t2, t11 * t47 - t148 * t67 - t1;];
tauc_reg = t5;
