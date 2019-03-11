% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:25
% EndTime: 2019-03-08 18:54:33
% DurationCPUTime: 3.33s
% Computational Cost: add. (4896->374), mult. (13252->537), div. (0->0), fcn. (11322->12), ass. (0->207)
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t160 = pkin(4) * t135 - pkin(10) * t138;
t110 = t160 * qJD(4);
t130 = sin(pkin(6));
t128 = sin(pkin(12));
t136 = sin(qJ(3));
t139 = cos(qJ(3));
t131 = cos(pkin(12));
t132 = cos(pkin(7));
t201 = t131 * t132;
t150 = t128 * t139 + t136 * t201;
t145 = t150 * t130;
t133 = cos(pkin(6));
t118 = qJD(1) * t133 + qJD(2);
t129 = sin(pkin(7));
t204 = t129 * t136;
t181 = t118 * t204;
t73 = qJD(1) * t145 + t181;
t258 = -t110 + t73;
t134 = sin(qJ(5));
t187 = qJD(5) * t134;
t172 = t135 * t187;
t137 = cos(qJ(5));
t185 = t137 * qJD(4);
t257 = -t138 * t185 + t172;
t191 = qJD(3) * t138;
t119 = -qJD(5) + t191;
t186 = qJD(5) * t137;
t70 = qJD(3) * pkin(9) + t73;
t194 = qJD(1) * t130;
t176 = t131 * t194;
t87 = t118 * t132 - t129 * t176;
t250 = -t135 * t70 + t138 * t87;
t203 = t129 * t139;
t205 = t128 * t136;
t247 = (t139 * t201 - t205) * t130;
t254 = qJD(1) * t247 + t118 * t203;
t67 = t254 * qJD(3);
t22 = qJD(4) * t250 + t138 * t67;
t224 = t135 * t87;
t47 = t138 * t70 + t224;
t44 = qJD(4) * pkin(10) + t47;
t55 = (t110 + t73) * qJD(3);
t113 = -pkin(4) * t138 - pkin(10) * t135 - pkin(3);
t72 = (t118 * t129 + t132 * t176) * t139 - t194 * t205;
t59 = t113 * qJD(3) - t72;
t167 = -t134 * t55 - t137 * t22 - t59 * t186 + t44 * t187;
t18 = -t134 * t44 + t137 * t59;
t256 = t119 * t18 - t167;
t255 = t133 * t203 + t247;
t19 = t134 * t59 + t137 * t44;
t7 = -qJD(5) * t19 - t134 * t22 + t137 * t55;
t253 = t19 * t119 - t7;
t189 = qJD(4) * t135;
t199 = t134 * t138;
t243 = pkin(9) * t134;
t252 = t258 * t137 - t189 * t243 - t72 * t199;
t197 = t137 * t138;
t251 = -t113 * t186 + t258 * t134 + t72 * t197;
t193 = qJD(3) * t135;
t105 = t134 * t193 - t185;
t211 = t105 * t119;
t182 = qJD(4) * qJD(5);
t84 = t257 * qJD(3) - t137 * t182;
t249 = -t84 + t211;
t190 = qJD(4) * t134;
t107 = t137 * t193 + t190;
t209 = t107 * t119;
t188 = qJD(4) * t138;
t148 = t134 * t188 + t135 * t186;
t85 = t148 * qJD(3) + t134 * t182;
t248 = t85 - t209;
t120 = pkin(9) * t197;
t90 = t134 * t113 + t120;
t246 = t107 ^ 2;
t245 = pkin(5) * t105;
t244 = pkin(5) * t134;
t23 = t47 * qJD(4) + t135 * t67;
t78 = t133 * t204 + t145;
t94 = -t129 * t130 * t131 + t132 * t133;
t57 = t135 * t78 - t94 * t138;
t242 = t23 * t57;
t95 = -t138 * t132 + t135 * t204;
t241 = t23 * t95;
t68 = t73 * qJD(3);
t240 = t68 * t255;
t239 = -qJ(6) - pkin(10);
t14 = -qJ(6) * t107 + t18;
t13 = -pkin(5) * t119 + t14;
t238 = t13 - t14;
t154 = pkin(5) * t135 - qJ(6) * t197;
t168 = qJD(5) * t239;
t109 = t160 * qJD(3);
t30 = t137 * t109 - t134 * t250;
t237 = t154 * qJD(3) + t134 * qJD(6) - t137 * t168 + t30;
t184 = t137 * qJD(6);
t31 = t134 * t109 + t137 * t250;
t236 = -t184 + t31 + (-qJ(6) * t191 - t168) * t134;
t235 = t135 * t184 - t154 * qJD(4) - (-t120 + (qJ(6) * t135 - t113) * t134) * qJD(5) + t252;
t234 = qJD(5) * t90 + t252;
t198 = t135 * t137;
t233 = -(-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t198 - (-qJD(6) * t135 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t138) * t134 + t251;
t232 = -(-t135 * t185 - t138 * t187) * pkin(9) + t251;
t231 = qJD(3) * pkin(3);
t230 = t105 * t72;
t229 = t107 * t72;
t16 = t85 * pkin(5) + t23;
t227 = t134 * t16;
t43 = -qJD(4) * pkin(4) - t250;
t226 = t134 * t43;
t225 = t135 * t72;
t223 = t137 * t16;
t222 = t137 * t43;
t220 = t139 * t68;
t15 = -qJ(6) * t105 + t19;
t219 = t15 * t119;
t217 = t23 * t134;
t216 = t23 * t135;
t215 = t23 * t137;
t214 = t84 * t134;
t213 = t85 * t137;
t210 = t107 * t105;
t208 = t107 * t134;
t207 = t119 * t134;
t206 = t119 * t137;
t141 = qJD(3) ^ 2;
t202 = t129 * t141;
t200 = t134 * t135;
t126 = t135 ^ 2;
t127 = t138 ^ 2;
t195 = t126 - t127;
t192 = qJD(3) * t136;
t183 = qJD(3) * qJD(4);
t179 = t136 * t202;
t177 = t135 * t141 * t138;
t175 = t129 * t192;
t174 = qJD(3) * t203;
t171 = t119 * t186;
t170 = t119 * t193;
t122 = t135 * t183;
t166 = -qJD(6) - t245;
t69 = -t72 - t231;
t165 = -qJD(3) * t69 - t67;
t164 = pkin(5) * t122;
t163 = t138 * t174;
t162 = t135 * t174;
t161 = t138 * t122;
t159 = -t134 * t19 - t137 * t18;
t58 = t135 * t94 + t138 * t78;
t29 = -t134 * t255 + t137 * t58;
t28 = -t134 * t58 - t137 * t255;
t157 = qJD(3) * t126 - t119 * t138;
t140 = qJD(4) ^ 2;
t156 = pkin(9) * t140;
t155 = qJD(4) * (t69 + t72 - t231);
t153 = t85 * qJ(6) + t167;
t96 = t132 * t135 + t138 * t204;
t81 = -t134 * t96 - t137 * t203;
t152 = t134 * t203 - t137 * t96;
t143 = t216 + t22 * t138 + (-t135 * t47 - t138 * t250) * qJD(4);
t142 = t84 * qJ(6) + t7;
t124 = -pkin(5) * t137 - pkin(4);
t115 = t239 * t137;
t114 = t239 * t134;
t111 = (pkin(9) + t244) * t135;
t104 = t137 * t113;
t102 = t105 ^ 2;
t89 = -pkin(9) * t199 + t104;
t88 = (-t119 - t191) * t189;
t86 = t148 * pkin(5) + pkin(9) * t188;
t83 = -qJ(6) * t200 + t90;
t80 = t96 * qJD(4) + t162;
t79 = -t95 * qJD(4) + t163;
t76 = -qJ(6) * t198 + t104 + (-pkin(5) - t243) * t138;
t75 = t78 * qJD(3);
t74 = t255 * qJD(3);
t71 = -t102 + t246;
t66 = -t209 - t85;
t65 = -t84 - t211;
t61 = -t171 + (t119 * t197 + (-t107 + t190) * t135) * qJD(3);
t60 = t119 * t187 + (-t119 * t199 + (t105 + t185) * t135) * qJD(3);
t51 = -t105 * t207 - t213;
t50 = -t107 * t206 - t214;
t49 = t148 * t105 + t85 * t200;
t48 = -t257 * t107 - t84 * t198;
t42 = t152 * qJD(5) - t134 * t79 + t137 * t175;
t41 = t81 * qJD(5) + t134 * t175 + t137 * t79;
t37 = t224 + (qJD(3) * t244 + t70) * t138;
t36 = t135 * t171 + t85 * t138 + (-t105 * t135 - t157 * t134) * qJD(4);
t35 = t119 * t172 + t84 * t138 + (t107 * t135 + t157 * t137) * qJD(4);
t32 = -t166 + t43;
t27 = -t57 * qJD(4) + t74 * t138;
t26 = t58 * qJD(4) + t74 * t135;
t20 = -t134 * t248 + t137 * t249;
t17 = (-t105 * t137 - t208) * t188 + (t214 - t213 + (t105 * t134 - t107 * t137) * qJD(5)) * t135;
t12 = t105 * t80 - t119 * t42 + t81 * t122 + t85 * t95;
t11 = t107 * t80 + t119 * t41 + t122 * t152 - t84 * t95;
t10 = t28 * qJD(5) + t75 * t134 + t27 * t137;
t9 = -t29 * qJD(5) - t27 * t134 + t75 * t137;
t8 = -t105 * t41 - t107 * t42 + t152 * t85 + t81 * t84;
t5 = -qJD(6) * t105 - t153;
t4 = -t107 * qJD(6) + t142 + t164;
t3 = t10 * t119 + t107 * t26 - t29 * t122 - t57 * t84;
t2 = t105 * t26 - t119 * t9 + t28 * t122 + t57 * t85;
t1 = -t10 * t105 - t107 * t9 + t28 * t84 - t29 * t85;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * qJD(3), -t74 * qJD(3), 0, t67 * t78 - t72 * t75 + t73 * t74 - t240, 0, 0, 0, 0, 0, 0, -t26 * qJD(4) + (-t138 * t75 - t189 * t255) * qJD(3), -t27 * qJD(4) + (t135 * t75 - t188 * t255) * qJD(3) (t135 * t26 + t138 * t27 + (-t135 * t58 + t138 * t57) * qJD(4)) * qJD(3), t22 * t58 - t250 * t26 + t27 * t47 + t69 * t75 - t240 + t242, 0, 0, 0, 0, 0, 0, t2, t3, t1, t10 * t19 - t167 * t29 + t18 * t9 + t26 * t43 + t28 * t7 + t242, 0, 0, 0, 0, 0, 0, t2, t3, t1, t10 * t15 + t13 * t9 + t16 * t57 + t26 * t32 + t28 * t4 + t29 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -t139 * t202, 0 (t136 * t67 - t220 + (-t136 * t72 + t139 * t73) * qJD(3)) * t129, 0, 0, 0, 0, 0, 0, -t138 * t179 + (-t80 - t162) * qJD(4), t135 * t179 + (-t79 - t163) * qJD(4) (t135 * t80 + t138 * t79 + (-t135 * t96 + t138 * t95) * qJD(4)) * qJD(3), t22 * t96 + t241 - t250 * t80 + t47 * t79 + (t69 * t192 - t220) * t129, 0, 0, 0, 0, 0, 0, t12, t11, t8, t152 * t167 + t18 * t42 + t19 * t41 + t43 * t80 + t7 * t81 + t241, 0, 0, 0, 0, 0, 0, t12, t11, t8, t13 * t42 + t15 * t41 - t152 * t5 + t16 * t95 + t32 * t80 + t4 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t150 * t194 - t181 + t73) * qJD(3) (-t254 + t72) * qJD(3), 0, 0, 0.2e1 * t161, -0.2e1 * t195 * t183, t140 * t138, -0.2e1 * t161, -t140 * t135, 0, t135 * t155 - t156 * t138, t156 * t135 + t138 * t155 (-t126 - t127) * t72 * qJD(3) + t143, -t68 * pkin(3) - t69 * t73 + (t135 * t250 - t138 * t47) * t72 + t143 * pkin(9), t48, t17, t35, t49, t36, t88, t234 * t119 + (-t7 + (pkin(9) * t105 + t226) * qJD(4)) * t138 + (t43 * t186 + pkin(9) * t85 - t230 + t217 + (qJD(3) * t89 + t18) * qJD(4)) * t135, -t232 * t119 + (-t167 + (pkin(9) * t107 + t222) * qJD(4)) * t138 + (-t43 * t187 - pkin(9) * t84 - t229 + t215 + (-qJD(3) * t90 - t19) * qJD(4)) * t135, t84 * t89 - t85 * t90 + t234 * t107 + t232 * t105 + t159 * t188 + (t134 * t167 - t137 * t7 + (t134 * t18 - t137 * t19) * qJD(5)) * t135, -t43 * t225 - t167 * t90 + t7 * t89 - t232 * t19 - t234 * t18 + (t188 * t43 + t216) * pkin(9), t48, t17, t35, t49, t36, t88, t105 * t86 + t111 * t85 + (t190 * t32 - t4) * t138 + t235 * t119 + (t32 * t186 - t230 + t227 + (qJD(3) * t76 + t13) * qJD(4)) * t135, t107 * t86 - t111 * t84 + (t185 * t32 + t5) * t138 - t233 * t119 + (-t32 * t187 - t229 + t223 + (-qJD(3) * t83 - t15) * qJD(4)) * t135, t76 * t84 - t83 * t85 + t235 * t107 + t233 * t105 + (-t13 * t137 - t134 * t15) * t188 + (-t134 * t5 - t137 * t4 + (t13 * t134 - t137 * t15) * qJD(5)) * t135, t16 * t111 + t4 * t76 + t5 * t83 + (t86 - t225) * t32 - t233 * t15 - t235 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, t195 * t141, 0, t177, 0, 0, t165 * t135, t165 * t138, 0, 0, t50, t20, t61, t51, t60, t170, -pkin(4) * t85 - t47 * t105 + t30 * t119 - t215 + (pkin(10) * t206 + t226) * qJD(5) + (-t135 * t18 + (-pkin(10) * t189 - t138 * t43) * t134) * qJD(3), pkin(4) * t84 - t47 * t107 - t31 * t119 + t217 + (-pkin(10) * t207 + t222) * qJD(5) + (-t43 * t197 + (-pkin(10) * t185 + t19) * t135) * qJD(3), t105 * t31 + t107 * t30 + ((qJD(5) * t107 - t85) * pkin(10) + t256) * t137 + ((qJD(5) * t105 - t84) * pkin(10) + t253) * t134, -t23 * pkin(4) - t18 * t30 - t19 * t31 - t43 * t47 + (qJD(5) * t159 - t7 * t134 - t137 * t167) * pkin(10), t50, t20, t61, t51, t60, t170, -t105 * t37 + t124 * t85 - t223 + t237 * t119 + (t32 + t245) * t187 + (-t32 * t199 + (qJD(4) * t114 - t13) * t135) * qJD(3), -t107 * t37 - t124 * t84 + t227 - t236 * t119 + (pkin(5) * t208 + t137 * t32) * qJD(5) + (-t32 * t197 + (qJD(4) * t115 + t15) * t135) * qJD(3), t114 * t84 + t115 * t85 + t237 * t107 + t236 * t105 + (t119 * t13 + t5) * t137 + (-t4 + t219) * t134, t4 * t114 - t5 * t115 + t16 * t124 + (pkin(5) * t187 - t37) * t32 - t236 * t15 - t237 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t71, t65, -t210, t66, t122, -t43 * t107 - t253, t105 * t43 - t256, 0, 0, t210, t71, t65, -t210, t66, t122, 0.2e1 * t164 - t219 + (t166 - t32) * t107 + t142, -t246 * pkin(5) - t14 * t119 + (qJD(6) + t32) * t105 + t153, t84 * pkin(5) - t238 * t105, t238 * t15 + (-t32 * t107 + t4) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t249, -t102 - t246, t15 * t105 + t13 * t107 + t16;];
tauc_reg  = t6;
