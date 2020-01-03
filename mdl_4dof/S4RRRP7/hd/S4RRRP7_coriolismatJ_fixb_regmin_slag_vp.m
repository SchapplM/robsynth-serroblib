% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:06
% EndTime: 2019-12-31 17:21:10
% DurationCPUTime: 1.79s
% Computational Cost: add. (1302->246), mult. (2911->350), div. (0->0), fcn. (2356->4), ass. (0->199)
t121 = sin(qJ(3));
t123 = cos(qJ(3));
t124 = cos(qJ(2));
t231 = pkin(5) * t121;
t122 = sin(qJ(2));
t228 = t122 * pkin(6);
t148 = -t124 * pkin(2) - t228;
t141 = -pkin(1) + t148;
t78 = t123 * t141;
t47 = -t78 + (pkin(3) + t231) * t124;
t179 = t124 * t231;
t55 = -t78 + t179;
t169 = -t55 / 0.2e1 + t47 / 0.2e1;
t226 = t123 * pkin(5);
t230 = t121 * pkin(2);
t46 = t121 * (-pkin(1) - t228) + (-qJ(4) + t226 - t230) * t124;
t203 = t123 * t124;
t178 = pkin(5) * t203;
t56 = t121 * t141 + t178;
t170 = -t46 / 0.2e1 + t56 / 0.2e1;
t246 = t170 * t121 + t169 * t123;
t117 = t121 ^ 2;
t119 = t123 ^ 2;
t101 = t119 - t117;
t187 = t122 * qJD(1);
t167 = t123 * t187;
t151 = t121 * t167;
t243 = t101 * qJD(2) - 0.2e1 * t151;
t207 = t121 * t122;
t108 = pkin(3) * t207;
t205 = t123 * qJ(4);
t59 = t108 + (pkin(5) - t205) * t122;
t238 = t59 / 0.2e1;
t208 = t121 * qJ(4);
t227 = t123 * pkin(3);
t146 = t208 + t227;
t71 = t146 * t122;
t237 = -t71 / 0.2e1;
t225 = t124 * pkin(6);
t229 = t122 * pkin(2);
t91 = -t225 + t229;
t236 = -t91 / 0.2e1;
t110 = pkin(5) * t207;
t235 = -t110 / 0.2e1;
t234 = -t122 / 0.2e1;
t233 = t122 / 0.2e1;
t232 = -t124 / 0.2e1;
t224 = t123 * t91;
t114 = t122 * qJ(4);
t206 = t122 * t123;
t83 = t121 * t91;
t49 = -pkin(5) * t206 + t114 + t83;
t115 = t122 * pkin(3);
t154 = -t110 - t224;
t50 = -t115 + t154;
t90 = pkin(3) * t121 - t205;
t60 = (pkin(5) + t90) * t124;
t3 = t46 * t49 + t47 * t50 + t59 * t60;
t223 = t3 * qJD(1);
t4 = -t46 * t55 + t47 * t56 + t59 * t71;
t222 = t4 * qJD(1);
t221 = t46 * t124;
t5 = (pkin(3) * t232 - t169) * t123 + (qJ(4) * t232 - t170) * t121;
t220 = t5 * qJD(1);
t219 = t55 * t124;
t218 = t56 * t124;
t217 = t59 * t121;
t216 = t59 * t123;
t7 = -t47 * t203 - t50 * t206 + (t122 * t49 + t221) * t121;
t215 = t7 * qJD(1);
t214 = t71 * t121;
t8 = -t56 * t206 + (t123 * t46 + (t47 - t55) * t121) * t122;
t213 = t8 * qJD(1);
t84 = -pkin(2) - t146;
t212 = t84 * t121;
t17 = t218 + (t214 + t216) * t122;
t211 = qJD(1) * t17;
t118 = t122 ^ 2;
t24 = -t118 * t231 - t219;
t210 = qJD(1) * t24;
t204 = t123 * t118;
t25 = -pkin(5) * t204 - t218;
t209 = qJD(1) * t25;
t144 = t122 * t60 + t124 * t59;
t13 = -t46 * t122 + t144 * t123 + t49 * t124;
t202 = t13 * qJD(1);
t14 = t144 * t121 - t47 * t122 + t50 * t124;
t201 = t14 * qJD(1);
t18 = t71 * t206 - t59 * t207 - t219;
t200 = t18 * qJD(1);
t19 = t55 * t122 + (-t154 - 0.2e1 * t110) * t124;
t199 = t19 * qJD(1);
t20 = t83 * t124 + (-t56 + t178) * t122;
t198 = t20 * qJD(1);
t21 = t59 * t206 + t221;
t197 = t21 * qJD(1);
t196 = t55 * qJD(3);
t120 = t124 ^ 2;
t102 = t120 - t118;
t81 = t102 * t121;
t195 = t81 * qJD(1);
t82 = t123 * t120 - t204;
t194 = t82 * qJD(1);
t193 = qJD(2) * t121;
t192 = qJD(2) * t123;
t191 = qJD(3) * t121;
t190 = qJD(3) * t123;
t189 = t102 * qJD(1);
t188 = t121 * qJD(4);
t186 = t122 * qJD(2);
t185 = t123 * qJD(4);
t184 = t124 * qJD(1);
t183 = t124 * qJD(2);
t182 = t124 * qJD(3);
t181 = t124 * qJD(4);
t180 = -t115 + t235;
t177 = pkin(1) * t187;
t176 = pkin(1) * t184;
t175 = pkin(6) * t191;
t174 = pkin(6) * t190;
t173 = t230 / 0.2e1;
t172 = -t225 / 0.2e1;
t171 = t225 / 0.2e1;
t168 = t84 * t233;
t166 = t121 * t182;
t165 = t121 * t190;
t105 = t121 * t192;
t164 = t121 * t185;
t163 = t122 * t183;
t162 = t122 * t190;
t161 = t122 * t188;
t160 = t122 * t184;
t159 = t123 * t186;
t158 = -t206 / 0.2e1;
t157 = t235 - t115 / 0.2e1;
t156 = t117 / 0.2e1 - t119 / 0.2e1;
t72 = t156 * t122;
t88 = t121 * qJD(1) * t204;
t155 = t72 * qJD(2) + t88;
t65 = (-0.1e1 / 0.2e1 + t156) * t122;
t153 = t65 * qJD(1) - t105;
t152 = t72 * qJD(1) - t105;
t107 = -qJD(3) + t184;
t150 = t121 * t159;
t149 = t90 * t233 + t238;
t147 = -t124 * t84 + t228;
t145 = t50 * t121 + t49 * t123;
t74 = -t83 / 0.2e1;
t93 = t121 * t171;
t10 = t93 - t214 / 0.2e1 - t216 / 0.2e1 - t114 + t74 + (t212 / 0.2e1 + (-t90 / 0.2e1 + pkin(5) / 0.2e1) * t123) * t122;
t37 = t90 * t121 + t84 * t123;
t143 = -t10 * qJD(1) + t37 * qJD(2);
t137 = t171 + t168;
t130 = t237 + t137;
t135 = t149 * t121;
t12 = t135 + (t236 + t130) * t123 + t180;
t38 = t90 * t123 - t212;
t142 = -t12 * qJD(1) + t38 * qJD(2);
t140 = t107 * t122;
t139 = t171 - t229 / 0.2e1;
t138 = -t49 * qJ(4) / 0.2e1 + t50 * pkin(3) / 0.2e1;
t125 = t246 * pkin(6) + t90 * t238 + t71 * t84 / 0.2e1;
t2 = t125 + t138;
t136 = -t84 * t90 * qJD(2) - t2 * qJD(1);
t34 = (t236 + t139) * t123;
t134 = pkin(2) * t193 - t34 * qJD(1);
t73 = t83 / 0.2e1;
t94 = t121 * t172;
t33 = t122 * t173 + t73 + t94;
t133 = pkin(2) * t192 - t33 * qJD(1);
t15 = t217 / 0.2e1 + (t236 + t137) * t123 + t157;
t132 = t15 * qJD(1) + t84 * t193;
t131 = t123 * t140;
t80 = t101 * t118;
t129 = t80 * qJD(1) + 0.2e1 * t150;
t127 = -t146 * qJD(3) + t185;
t87 = t119 * t118 + t120;
t126 = t87 * qJD(1) + t150 - t182;
t111 = t186 / 0.2e1;
t106 = t123 * t184;
t89 = t123 * t161;
t86 = t107 * qJ(4);
t79 = -t106 + t190;
t77 = (t184 - qJD(3) / 0.2e1) * t122;
t75 = t224 / 0.2e1;
t70 = t117 * qJD(2) + t151;
t69 = t72 * qJD(3);
t66 = t117 * t234 + t119 * t233 + t234;
t62 = (-t121 * t187 + t192) * t124;
t61 = t121 * t140;
t54 = t56 * qJD(3);
t23 = t139 * t123 + t110 + t75;
t22 = t94 + t74 + (t173 + t226) * t122;
t16 = t123 * t172 + t84 * t158 - t217 / 0.2e1 - t224 / 0.2e1 + t157;
t11 = t130 * t123 + t135 - t180 + t75;
t9 = t93 + t114 + pkin(5) * t158 + t73 - t149 * t123 + (t168 + t237) * t121;
t6 = (-t208 / 0.2e1 - t227 / 0.2e1) * t124 + t246;
t1 = t125 - t138;
t26 = [0, 0, 0, t163, t102 * qJD(2), 0, 0, 0, -pkin(1) * t186, -pkin(1) * t183, -t118 * t165 + t119 * t163, -t80 * qJD(3) - 0.2e1 * t124 * t150, -t82 * qJD(2) + t122 * t166, t81 * qJD(2) + t124 * t162, -t163, -qJD(2) * t19 - qJD(3) * t25, qJD(2) * t20 + qJD(3) * t24, t14 * qJD(2) + t17 * qJD(3) - t118 * t164, -t7 * qJD(2) - t8 * qJD(3) + t124 * t161, -t13 * qJD(2) - t18 * qJD(3) + t87 * qJD(4), qJD(2) * t3 + qJD(3) * t4 - qJD(4) * t21; 0, 0, 0, t160, t189, t183, -t186, 0, -pkin(5) * t183 - t177, pkin(5) * t186 - t176, -t69 + (t119 * t187 + t105) * t124, -0.2e1 * t121 * t162 + t243 * t124, t121 * t186 - t194, t159 + t195, -t77, -t199 + (t121 * t148 - t178) * qJD(2) + t23 * qJD(3), t198 + (t123 * t148 + t179) * qJD(2) + t22 * qJD(3), t201 + (-t121 * t147 - t60 * t123) * qJD(2) + t11 * qJD(3) + t66 * qJD(4), qJD(2) * t145 + t6 * qJD(3) - t215, -t202 + (-t60 * t121 + t123 * t147) * qJD(2) + t9 * qJD(3) + t89, t223 + (pkin(6) * t145 + t60 * t84) * qJD(2) + t1 * qJD(3) + t16 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t129, t61, t131, t111, qJD(2) * t23 - t209 - t54, qJD(2) * t22 + t196 + t210, qJD(2) * t11 + t211 - t54, -t213 + t6 * qJD(2) + (-t122 * t205 + t108) * qJD(3) - t161, t9 * qJD(2) - t181 - t196 - t200, t222 + t1 * qJD(2) + (-pkin(3) * t56 - qJ(4) * t55) * qJD(3) + t46 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * qJD(2) - t88, t61, t126, qJD(2) * t16 + qJD(3) * t46 - t197; 0, 0, 0, -t160, -t189, 0, 0, 0, t177, t176, -t119 * t160 - t69, 0.2e1 * t121 * t131, -t123 * t182 + t194, t166 - t195, t77, qJD(3) * t34 + t199, qJD(3) * t33 - t198, t12 * qJD(3) - t65 * qJD(4) - t201, -t5 * qJD(3) - t123 * t181 + t215, t10 * qJD(3) + t202 + t89, qJD(3) * t2 - qJD(4) * t15 - t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t101 * qJD(3), 0, 0, 0, -pkin(2) * t191, -pkin(2) * t190, -t38 * qJD(3) + t164, 0, -t37 * qJD(3) + t117 * qJD(4), (qJD(3) * t90 - t188) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t243, t79, t107 * t121, -t187 / 0.2e1, -t134 - t174, -t133 + t175, -t142 - t174, t127 - t220, -t143 - t175, pkin(6) * t127 - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t79, t70, -t132 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t129, t62, (-t167 - t193) * t124, t111, -qJD(2) * t34 + t209, -qJD(2) * t33 - t210, -qJD(2) * t12 - t211, qJD(2) * t5 + t213, -t10 * qJD(2) - t181 + t200, -qJ(4) * t181 - t2 * qJD(2) - t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t243, t106, -t121 * t184, t187 / 0.2e1, t134, t133, t142, t220, t143, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * qJD(2) + t88, t62, -t126, qJ(4) * t182 + t15 * qJD(2) + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t106, -t70, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t26;
