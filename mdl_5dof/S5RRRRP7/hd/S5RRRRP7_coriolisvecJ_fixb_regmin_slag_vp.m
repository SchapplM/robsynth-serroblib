% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:52
% DurationCPUTime: 2.37s
% Computational Cost: add. (4087->319), mult. (10027->425), div. (0->0), fcn. (6866->6), ass. (0->177)
t128 = sin(qJ(2));
t127 = sin(qJ(3));
t130 = cos(qJ(2));
t191 = t127 * t130;
t221 = cos(qJ(3));
t104 = t128 * t221 + t191;
t179 = qJD(2) + qJD(3);
t75 = t179 * t104;
t233 = t75 * qJD(1);
t126 = sin(qJ(4));
t182 = qJD(4) * t126;
t168 = qJD(1) * t221;
t185 = qJD(1) * t128;
t94 = t127 * t185 - t130 * t168;
t201 = t126 * t94;
t232 = -t182 - t201;
t180 = qJD(1) * qJD(2);
t231 = -0.2e1 * t180;
t129 = cos(qJ(4));
t181 = qJD(4) * t129;
t144 = -t127 * t128 + t130 * t221;
t140 = t144 * qJD(3);
t29 = t233 * pkin(3) + (-pkin(8) * t140 + (t128 * pkin(2) - pkin(8) * t144) * qJD(2)) * qJD(1);
t224 = -pkin(7) - pkin(6);
t175 = qJD(2) * t224;
t159 = qJD(1) * t175;
t100 = t128 * t159;
t101 = t130 * t159;
t113 = t224 * t130;
t108 = qJD(1) * t113;
t167 = t221 * qJD(3);
t184 = qJD(3) * t127;
t112 = t224 * t128;
t106 = qJD(1) * t112;
t206 = qJD(2) * pkin(2);
t99 = t106 + t206;
t36 = t100 * t221 + t127 * t101 + t108 * t184 + t167 * t99;
t122 = -pkin(2) * t130 - pkin(1);
t111 = t122 * qJD(1);
t96 = -qJD(1) * t191 - t128 * t168;
t57 = t94 * pkin(3) + t96 * pkin(8) + t111;
t98 = t221 * t108;
t71 = t127 * t99 - t98;
t60 = pkin(8) * t179 + t71;
t146 = t126 * t29 + t129 * t36 + t181 * t57 - t182 * t60;
t193 = t233 * qJ(5);
t91 = qJD(4) + t94;
t2 = qJD(5) * t91 + t146 + t193;
t165 = t126 * t36 - t129 * t29 + t181 * t60 + t182 * t57;
t223 = t233 * pkin(4);
t5 = t165 - t223;
t230 = t5 * t126 + t2 * t129;
t199 = t129 * t233;
t147 = t182 * t91 - t199;
t229 = pkin(8) * t147;
t192 = qJ(5) * t129;
t228 = pkin(4) * t232 + qJ(5) * t181 + t126 * qJD(5) + t94 * t192;
t72 = t106 * t127 - t98;
t158 = pkin(2) * t184 - t72;
t227 = t221 * t112 + t113 * t127;
t74 = qJD(2) * t144 + t140;
t134 = t74 * qJD(1);
t143 = -t126 * t179 + t129 * t96;
t42 = -qJD(4) * t143 + t126 * t134;
t226 = t143 ^ 2;
t225 = t91 ^ 2;
t222 = t96 * pkin(4);
t97 = t127 * t108;
t70 = t221 * t99 + t97;
t59 = -pkin(3) * t179 - t70;
t162 = t129 * t179;
t77 = -t126 * t96 - t162;
t22 = pkin(4) * t77 + qJ(5) * t143 + t59;
t220 = t22 * t143;
t219 = t59 * t94;
t37 = t100 * t127 - t101 * t221 - t108 * t167 + t184 * t99;
t41 = -qJD(4) * t162 - t129 * t134 - t182 * t96;
t7 = pkin(4) * t42 + qJ(5) * t41 + qJD(5) * t143 + t37;
t218 = t7 * t129;
t217 = t77 * t91;
t216 = t143 * t77;
t215 = t143 * t91;
t214 = t91 * t96;
t213 = t96 * t94;
t212 = t71 + t228;
t211 = -t158 + t228;
t67 = -pkin(3) * t96 + pkin(8) * t94;
t210 = t126 * t67 + t129 * t70;
t58 = pkin(2) * t185 + t67;
t73 = t106 * t221 + t97;
t209 = t126 * t58 + t129 * t73;
t69 = -pkin(3) * t144 - pkin(8) * t104 + t122;
t82 = t112 * t127 - t113 * t221;
t208 = t126 * t69 + t129 * t82;
t207 = pkin(2) * qJD(3);
t120 = pkin(2) * t127 + pkin(8);
t205 = t120 * t233;
t204 = t126 * t233;
t203 = t126 * t74;
t202 = t126 * t77;
t200 = t129 * t42;
t198 = t129 * t74;
t197 = t129 * t143;
t164 = t129 * t91;
t196 = t129 * t94;
t195 = t37 * t129;
t194 = t41 * t126;
t132 = qJD(1) ^ 2;
t190 = t130 * t132;
t131 = qJD(2) ^ 2;
t189 = t131 * t128;
t188 = t131 * t130;
t25 = -t126 * t60 + t129 * t57;
t187 = qJD(5) - t25;
t186 = t128 ^ 2 - t130 ^ 2;
t183 = qJD(4) * t120;
t178 = t221 * pkin(2);
t177 = t128 * t206;
t19 = t22 * t182;
t53 = t59 * t182;
t174 = t126 * t221;
t173 = t129 * t221;
t26 = t126 * t57 + t129 * t60;
t18 = qJ(5) * t91 + t26;
t171 = t18 * t96 - t196 * t22;
t17 = -pkin(4) * t91 + t187;
t170 = -t17 * t96 + t19;
t169 = t25 * t96 + t53;
t166 = t128 * t180;
t163 = pkin(1) * t231;
t161 = pkin(2) * t167;
t160 = t126 * t37 + t181 * t59 - t26 * t96;
t156 = -t183 * t91 - t7;
t155 = pkin(4) * t129 + qJ(5) * t126;
t154 = pkin(4) * t126 - t192;
t153 = -t126 * t18 + t129 * t17;
t152 = -t126 * t70 + t129 * t67;
t110 = -pkin(3) - t155;
t151 = t26 * t91 - t165;
t150 = t7 * t126 + t181 * t22;
t149 = t111 * t96 - t37;
t148 = -t181 * t91 - t204;
t40 = pkin(3) * t75 - pkin(8) * t74 + t177;
t107 = t128 * t175;
t109 = t130 * t175;
t45 = qJD(3) * t227 + t107 * t221 + t127 * t109;
t145 = t126 * t40 + t129 * t45 + t181 * t69 - t182 * t82;
t142 = t148 * pkin(8);
t141 = -t161 * t91 - t205;
t139 = t230 + t232 * t18 + (t181 + t196) * t17;
t138 = qJD(4) * t153 + t230;
t137 = -t194 - t200 + (-t197 + t202) * qJD(4);
t136 = t111 * t94 - t36;
t46 = qJD(3) * t82 + t127 * t107 - t109 * t221;
t121 = -t178 - pkin(3);
t102 = -t178 + t110;
t89 = t96 * qJ(5);
t50 = -t94 ^ 2 + t96 ^ 2;
t49 = -t179 * t96 - t233;
t48 = t179 * t94 + t134;
t47 = -pkin(4) * t143 + qJ(5) * t77;
t44 = t104 * t154 - t227;
t31 = pkin(4) * t144 + t126 * t82 - t129 * t69;
t30 = -qJ(5) * t144 + t208;
t24 = -t152 + t222;
t23 = -t89 + t210;
t21 = t126 * t73 - t129 * t58 + t222;
t20 = -t89 + t209;
t16 = -t41 + t217;
t12 = -t143 * t96 + t164 * t91 + t204;
t11 = -t126 * t225 - t77 * t96 + t199;
t10 = -t143 * t164 - t194;
t9 = t154 * t74 + (qJD(4) * t155 - qJD(5) * t129) * t104 + t46;
t8 = -t75 * pkin(4) + qJD(4) * t208 + t126 * t45 - t129 * t40;
t6 = qJ(5) * t75 - qJD(5) * t144 + t145;
t3 = (-t41 - t217) * t129 + (-t42 + t215) * t126;
t1 = [0, 0, 0, 0.2e1 * t130 * t166, t186 * t231, t188, -t189, 0, -pkin(6) * t188 + t128 * t163, pkin(6) * t189 + t130 * t163, t104 * t134 - t96 * t74, -t104 * t233 + t134 * t144 - t74 * t94 + t96 * t75, t74 * t179, -t75 * t179, 0, t122 * t233 + t111 * t75 - t46 * t179 + (-qJD(1) * t144 + t94) * t177, pkin(2) * t104 * t166 + t111 * t74 + t122 * t134 - t177 * t96 - t179 * t45, -t74 * t197 + (-t41 * t129 + t143 * t182) * t104, (t126 * t143 - t129 * t77) * t74 + (t194 - t200 + (t197 + t202) * qJD(4)) * t104, -t104 * t147 - t143 * t75 + t144 * t41 + t164 * t74, t104 * t148 + t144 * t42 - t203 * t91 - t77 * t75, -t144 * t233 + t75 * t91, t165 * t144 + t25 * t75 + t46 * t77 - t227 * t42 + ((-qJD(4) * t82 + t40) * t91 + t69 * t233 + t59 * qJD(4) * t104) * t129 + ((-qJD(4) * t69 - t45) * t91 - t82 * t233 + t37 * t104 + t59 * t74) * t126, -t145 * t91 - t208 * t233 + t146 * t144 - t26 * t75 - t46 * t143 + t227 * t41 + t59 * t198 + (-t53 + t195) * t104, t104 * t150 + t144 * t5 - t17 * t75 + t203 * t22 - t233 * t31 + t44 * t42 + t9 * t77 - t8 * t91, -t30 * t42 - t31 * t41 - t6 * t77 - t8 * t143 + t153 * t74 + (-t126 * t2 + t129 * t5 + (-t126 * t17 - t129 * t18) * qJD(4)) * t104, -t22 * t198 - t2 * t144 + t18 * t75 + t30 * t233 + t44 * t41 + t6 * t91 + t9 * t143 + (t19 - t218) * t104, t17 * t8 + t18 * t6 + t2 * t30 + t22 * t9 + t31 * t5 + t44 * t7; 0, 0, 0, -t128 * t190, t186 * t132, 0, 0, 0, t132 * pkin(1) * t128, pkin(1) * t190, -t213, t50, t48, t49, 0, t72 * t179 + (-t179 * t184 - t185 * t94) * pkin(2) + t149, t73 * t179 + (-t167 * t179 + t185 * t96) * pkin(2) + t136, t10, t3, t12, t11, t214, t121 * t42 + t158 * t77 + (-t37 + (-t58 - t183) * t91) * t129 + (-t205 + t219 + (-t161 + t73) * t91) * t126 + t169, -t121 * t41 + (t120 * t182 + t209) * t91 - t158 * t143 + (t141 + t219) * t129 + t160, t102 * t42 + t21 * t91 - t211 * t77 + t156 * t129 + (t22 * t94 + t141) * t126 + t170, t20 * t77 + t21 * t143 + (-t143 * t174 - t173 * t77) * t207 + t137 * t120 + t139, t102 * t41 - t20 * t91 - t211 * t143 + t156 * t126 + (-qJD(4) * t22 - t141) * t129 + t171, t7 * t102 - t17 * t21 - t18 * t20 - t211 * t22 + (t17 * t174 + t173 * t18) * t207 + t138 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t50, t48, t49, 0, t179 * t71 + t149, t179 * t70 + t136, t10, t3, t12, t11, t214, -pkin(3) * t42 - t152 * t91 + t201 * t59 - t71 * t77 + t142 + t169 - t195, pkin(3) * t41 + t143 * t71 + t196 * t59 + t210 * t91 + t160 + t229, t110 * t42 + t201 * t22 - t212 * t77 + t24 * t91 + t142 + t170 - t218, pkin(8) * t137 + t143 * t24 + t23 * t77 + t139, t110 * t41 - t143 * t212 - t23 * t91 - t150 + t171 - t229, pkin(8) * t138 + t7 * t110 - t17 * t24 - t18 * t23 - t212 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216, -t77 ^ 2 + t226, t16, -t42 - t215, t233, t143 * t59 + t151, t25 * t91 + t59 * t77 - t146, -t47 * t77 + t151 + t220 + 0.2e1 * t223, pkin(4) * t41 - t42 * qJ(5) - (t18 - t26) * t143 + (t17 - t187) * t77, 0.2e1 * t193 - t22 * t77 - t47 * t143 + (0.2e1 * qJD(5) - t25) * t91 + t146, -t5 * pkin(4) + t2 * qJ(5) - t17 * t26 + t18 * t187 - t22 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233 - t216, t16, -t225 - t226, -t18 * t91 - t220 + t5;];
tauc_reg = t1;
