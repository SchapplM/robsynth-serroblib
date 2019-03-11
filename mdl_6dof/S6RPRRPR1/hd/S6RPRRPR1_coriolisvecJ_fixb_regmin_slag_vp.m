% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:45
% EndTime: 2019-03-09 04:58:52
% DurationCPUTime: 2.26s
% Computational Cost: add. (4864->279), mult. (11770->383), div. (0->0), fcn. (8622->10), ass. (0->178)
t142 = cos(qJ(4));
t143 = cos(qJ(3));
t201 = qJD(1) * t143;
t190 = t142 * t201;
t139 = sin(qJ(4));
t140 = sin(qJ(3));
t202 = qJD(1) * t140;
t191 = t139 * t202;
t102 = -t190 + t191;
t104 = -t139 * t201 - t142 * t202;
t135 = sin(pkin(11));
t212 = cos(pkin(11));
t176 = -t102 * t212 + t104 * t135;
t236 = qJD(6) - t176;
t245 = qJD(6) - t236;
t124 = sin(pkin(10)) * pkin(1) + pkin(7);
t225 = pkin(8) + t124;
t132 = qJD(3) + qJD(4);
t138 = sin(qJ(6));
t141 = cos(qJ(6));
t158 = -t102 * t135 - t104 * t212;
t58 = -t132 * t141 + t138 * t158;
t244 = t236 * t58;
t175 = t141 * t236;
t109 = t139 * t143 + t140 * t142;
t237 = qJD(1) * t109;
t148 = t132 * t237;
t196 = qJD(1) * qJD(3);
t187 = t143 * t196;
t68 = qJD(4) * t190 - t132 * t191 + t142 * t187;
t44 = t135 * t68 + t148 * t212;
t217 = t138 * t44;
t243 = -t175 * t236 - t217;
t129 = pkin(3) * t202;
t233 = t104 * pkin(4);
t37 = pkin(5) * t158 - pkin(9) * t176 - t233;
t128 = pkin(3) * t142 + pkin(4);
t181 = t212 * t139;
t97 = pkin(3) * t181 + t128 * t135;
t93 = pkin(9) + t97;
t242 = (qJD(6) * t93 + t129 + t37) * t236;
t180 = t225 * qJD(1);
t85 = t140 * qJD(2) + t143 * t180;
t80 = t139 * t85;
t219 = qJD(3) * pkin(3);
t84 = t143 * qJD(2) - t140 * t180;
t83 = t84 + t219;
t185 = t142 * t83 - t80;
t98 = t104 * qJ(5);
t49 = t185 + t98;
t199 = qJD(6) * t138;
t240 = -t141 * t44 + t199 * t236;
t78 = t84 * qJD(3);
t239 = (qJD(4) * t83 + t78) * t142;
t108 = t139 * t140 - t142 * t143;
t73 = -t108 * t135 + t109 * t212;
t226 = t73 * t44;
t75 = t132 * t108;
t156 = t109 * qJD(4);
t76 = qJD(3) * t109 + t156;
t53 = -t135 * t76 - t212 * t75;
t169 = t236 * t53 + t226;
t194 = t73 * t199;
t235 = -t141 * t169 + t194 * t236;
t182 = qJD(3) * t225;
t100 = t143 * t182;
t106 = t225 * t140;
t107 = t225 * t143;
t164 = t106 * t139 - t107 * t142;
t99 = t140 * t182;
t152 = qJD(4) * t164 - t100 * t142 + t139 * t99;
t149 = t75 * qJ(5) - t109 * qJD(5) + t152;
t200 = qJD(4) * t139;
t207 = t142 * t106;
t159 = -qJD(4) * t207 - t100 * t139 - t107 * t200 - t142 * t99;
t31 = -qJ(5) * t76 - qJD(5) * t108 + t159;
t11 = t135 * t149 + t212 * t31;
t82 = t142 * t85;
t166 = -t139 * t83 - t82;
t79 = t85 * qJD(3);
t167 = -t139 * t78 - t142 * t79;
t150 = qJD(4) * t166 + t167;
t147 = -t68 * qJ(5) + t104 * qJD(5) + t150;
t184 = -t139 * t79 - t200 * t85;
t15 = -qJ(5) * t148 - t102 * qJD(5) + t184 + t239;
t3 = t135 * t15 - t147 * t212;
t154 = -qJ(5) * t109 - t107 * t139 - t207;
t57 = -qJ(5) * t108 - t164;
t35 = t135 * t154 + t212 * t57;
t172 = t3 * t73 - t35 * t44;
t211 = t102 * qJ(5);
t50 = -t166 - t211;
t218 = t135 * t50;
t43 = pkin(4) * t132 + t49;
t22 = t212 * t43 - t218;
t19 = -pkin(5) * t132 - t22;
t126 = -cos(pkin(10)) * pkin(1) - pkin(2);
t110 = -pkin(3) * t143 + t126;
t105 = t110 * qJD(1);
t71 = t102 * pkin(4) + qJD(5) + t105;
t33 = -pkin(5) * t176 - pkin(9) * t158 + t71;
t155 = pkin(4) * t108 + t110;
t72 = t108 * t212 + t109 * t135;
t38 = pkin(5) * t72 - pkin(9) * t73 + t155;
t4 = t135 * t147 + t15 * t212;
t234 = -(qJD(6) * t38 + t11) * t236 - (qJD(6) * t33 + t4) * t72 + t19 * t53 + t172;
t232 = t19 * t176;
t231 = t19 * t73;
t230 = t38 * t44;
t60 = t132 * t138 + t141 * t158;
t229 = t60 * t158;
t228 = t236 * t158;
t227 = t158 * t58;
t198 = qJD(6) * t141;
t45 = -t135 * t148 + t212 * t68;
t28 = t132 * t198 + t141 * t45 - t158 * t199;
t52 = -t135 * t75 + t212 * t76;
t224 = t28 * t72 + t52 * t60;
t46 = t212 * t50;
t23 = t135 * t43 + t46;
t223 = t142 * t84 - t80;
t183 = -t139 * t84 - t82;
t161 = t183 + t211;
t220 = pkin(3) * qJD(4);
t51 = t98 + t223;
t222 = -t135 * t51 + t212 * t161 + (t135 * t142 + t181) * t220;
t208 = t135 * t139;
t221 = -t135 * t161 - t212 * t51 + (t142 * t212 - t208) * t220;
t216 = t138 * t45;
t215 = t138 * t176;
t214 = t28 * t138;
t213 = t75 * t132;
t210 = t104 * t102;
t209 = t105 * t104;
t144 = qJD(3) ^ 2;
t206 = t144 * t140;
t205 = t144 * t143;
t203 = t140 ^ 2 - t143 ^ 2;
t112 = qJD(1) * t126;
t130 = t140 * t219;
t193 = t236 * t198;
t20 = pkin(9) * t132 + t23;
t168 = t138 * t20 - t141 * t33;
t192 = t158 * t168 + t19 * t199;
t189 = -pkin(3) * t132 - t83;
t188 = pkin(4) * t76 + t130;
t8 = t138 * t33 + t141 * t20;
t174 = t138 * t3 + t158 * t8 + t19 * t198;
t171 = t158 * t23 + t176 * t22;
t29 = qJD(6) * t60 + t216;
t170 = -t72 * t29 - t52 * t58;
t163 = t215 * t236 - t240;
t162 = 0.2e1 * qJD(3) * t112;
t160 = t105 * t102 - t184;
t96 = -pkin(3) * t208 + t128 * t212;
t153 = -t221 * t236 - t93 * t44 - t232;
t151 = -t138 * t169 - t193 * t73;
t146 = pkin(4) * t148 + qJD(3) * t129;
t145 = qJD(1) ^ 2;
t125 = -pkin(4) * t212 - pkin(5);
t123 = pkin(4) * t135 + pkin(9);
t92 = -pkin(5) - t96;
t69 = t76 * t132;
t61 = -t102 ^ 2 + t104 ^ 2;
t55 = -t104 * t132 - t148;
t54 = t102 * t132 + t68;
t34 = t135 * t57 - t154 * t212;
t25 = t212 * t49 - t218;
t24 = t135 * t49 + t46;
t16 = pkin(5) * t52 - pkin(9) * t53 + t188;
t13 = t44 * pkin(5) - t45 * pkin(9) + t146;
t12 = t141 * t13;
t10 = t135 * t31 - t149 * t212;
t9 = t175 * t60 + t214;
t6 = -t229 - t243;
t5 = t163 + t227;
t1 = (t28 - t244) * t141 + (-t236 * t60 - t29) * t138;
t2 = [0, 0, 0, 0, 0.2e1 * t140 * t187, -0.2e1 * t203 * t196, t205, -t206, 0, -t124 * t205 + t140 * t162, t124 * t206 + t143 * t162, t104 * t75 + t109 * t68, t75 * t102 + t104 * t76 - t68 * t108 - t109 * t148, -t213, -t69, 0, t102 * t130 + t105 * t76 + t152 * t132 + (t110 * t156 + (t140 * pkin(3) * t108 + t109 * t110) * qJD(3)) * qJD(1), t110 * t68 - t105 * t75 - t159 * t132 + (-t104 + t237) * t130, t10 * t158 + t11 * t176 - t22 * t53 - t23 * t52 + t34 * t45 - t4 * t72 + t172, -t22 * t10 + t23 * t11 + t146 * t155 + t188 * t71 + t3 * t34 + t4 * t35, -t60 * t194 + (t28 * t73 + t53 * t60) * t141 (-t138 * t60 - t141 * t58) * t53 + (-t214 - t141 * t29 + (t138 * t58 - t141 * t60) * qJD(6)) * t73, t224 - t235, t151 + t170, t236 * t52 + t44 * t72, t10 * t58 + t12 * t72 + t34 * t29 - t168 * t52 + (t16 * t236 + t230 + (-t20 * t72 - t236 * t35 + t231) * qJD(6)) * t141 + t234 * t138, t10 * t60 + t34 * t28 - t8 * t52 + (-(-qJD(6) * t35 + t16) * t236 - t230 - (-qJD(6) * t20 + t13) * t72 - qJD(6) * t231) * t138 + t234 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206, -t205, 0, 0, 0, 0, 0, -t69, t213, t158 * t52 + t176 * t53 + t45 * t72 - t226, -t22 * t52 + t23 * t53 + t3 * t72 + t4 * t73, 0, 0, 0, 0, 0, t151 - t170, t224 + t235; 0, 0, 0, 0, -t140 * t145 * t143, t203 * t145, 0, 0, 0, -t112 * t202, -t112 * t201, -t210, t61, t54, t55, 0, -t102 * t129 + t209 - t183 * t132 + (t139 * t189 - t82) * qJD(4) + t167, t104 * t129 + t223 * t132 + (qJD(4) * t189 - t78) * t142 + t160, t158 * t222 + t176 * t221 - t97 * t44 - t96 * t45 + t171, t4 * t97 - t3 * t96 - t71 * (t129 - t233) + t221 * t23 - t222 * t22, t9, t1, t6, t5, -t228, t92 * t29 + t222 * t58 + (-t3 - t242) * t141 + t153 * t138 + t192, t138 * t242 + t141 * t153 + t222 * t60 + t92 * t28 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, t61, t54, t55, 0, -t132 * t166 + t150 + t209, t132 * t185 + t160 - t239, -t24 * t158 - t25 * t176 + (-t135 * t44 - t212 * t45) * pkin(4) + t171, t22 * t24 - t23 * t25 + (t104 * t71 + t135 * t4 - t212 * t3) * pkin(4), t9, t1, t6, t5, -t228, t125 * t29 - t3 * t141 - (-t138 * t25 + t141 * t37) * t236 - t24 * t58 - t19 * t215 + (-t193 - t217) * t123 + t192, t125 * t28 + (t138 * t37 + t141 * t25) * t236 - t24 * t60 - t141 * t232 + t240 * t123 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158 ^ 2 - t176 ^ 2, t158 * t22 - t176 * t23 + t146, 0, 0, 0, 0, 0, t163 - t227, -t229 + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, t28 + t244, -t245 * t60 - t216, t44, -t138 * t4 - t19 * t60 - t245 * t8 + t12, -t138 * t13 - t141 * t4 + t168 * t245 + t19 * t58;];
tauc_reg  = t2;
