% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:36:58
% EndTime: 2019-03-09 04:37:07
% DurationCPUTime: 3.02s
% Computational Cost: add. (2982->370), mult. (6688->456), div. (0->0), fcn. (3936->6), ass. (0->182)
t132 = cos(qJ(3));
t205 = qJD(1) * t132;
t112 = -qJD(4) + t205;
t130 = sin(qJ(3));
t124 = t130 ^ 2;
t131 = cos(qJ(4));
t197 = qJD(1) * qJD(3);
t180 = t131 * t197;
t129 = sin(qJ(4));
t200 = qJD(4) * t129;
t182 = t130 * t200;
t198 = t131 * qJD(3);
t184 = t132 * t198;
t146 = (t182 - t184) * t112 + t124 * t180;
t202 = qJD(3) * t130;
t196 = qJD(3) * qJD(4);
t47 = qJD(1) * t182 - t131 * t196 - t132 * t180;
t203 = qJD(3) * t129;
t206 = qJD(1) * t130;
t87 = t131 * t206 + t203;
t252 = -t132 * t47 - t87 * t202;
t260 = t146 + t252;
t117 = t130 * t197;
t240 = pkin(4) + qJ(6);
t259 = t240 * t117;
t199 = qJD(4) * t131;
t181 = t130 * t199;
t137 = -(qJD(1) * t124 - t112 * t132) * t203 + t112 * t181;
t201 = qJD(3) * t132;
t258 = t129 * t201 + t181;
t48 = qJD(1) * t258 + t129 * t196;
t85 = t129 * t206 - t198;
t173 = t132 * t48 - t85 * t202;
t135 = t137 - t173;
t121 = t130 * qJD(2);
t115 = sin(pkin(9)) * pkin(1) + pkin(7);
t97 = t115 * qJD(1);
t65 = t132 * t97 + t121;
t55 = qJD(3) * pkin(8) + t65;
t116 = -cos(pkin(9)) * pkin(1) - pkin(2);
t78 = -pkin(3) * t132 - pkin(8) * t130 + t116;
t58 = t78 * qJD(1);
t23 = t129 * t55 - t131 * t58;
t153 = pkin(5) * t87 + t23;
t209 = qJD(5) + t153;
t24 = t129 * t58 + t131 * t55;
t18 = qJ(5) * t112 - t24;
t166 = pkin(4) * t117;
t254 = qJD(2) * t132 - t130 * t97;
t56 = t254 * qJD(3);
t162 = pkin(3) * t130 - pkin(8) * t132;
t94 = t162 * qJD(3);
t77 = qJD(1) * t94;
t171 = t129 * t56 - t131 * t77 + t55 * t199 + t58 * t200;
t5 = -t166 + t171;
t256 = -t112 * t18 + t5;
t110 = t112 ^ 2;
t82 = t87 ^ 2;
t255 = -t82 - t110;
t102 = qJD(5) * t112;
t111 = qJ(5) * t117;
t253 = t111 - t102;
t251 = qJD(5) * t129 + t65 + (t129 * t205 - t200) * pkin(4);
t249 = t85 ^ 2;
t248 = pkin(5) + pkin(8);
t247 = pkin(5) * t85;
t9 = t240 * t112 + t209;
t246 = t112 * t9;
t57 = qJD(3) * t121 + t97 * t201;
t144 = qJ(5) * t47 - qJD(5) * t87 + t57;
t8 = pkin(4) * t48 + t144;
t245 = t129 * t8;
t244 = t131 * t8;
t54 = -qJD(3) * pkin(3) - t254;
t140 = -qJ(5) * t87 + t54;
t14 = t240 * t85 + t140;
t243 = t14 * t87;
t22 = pkin(4) * t85 + t140;
t242 = t22 * t87;
t241 = t87 * t85;
t101 = t248 * t131;
t91 = t162 * qJD(1);
t175 = -t129 * t254 + t131 * t91;
t213 = t131 * t132;
t193 = pkin(5) * t213;
t239 = -(-t240 * t130 + t193) * qJD(1) + t175 + qJD(4) * t101;
t215 = t129 * t132;
t194 = pkin(5) * t215;
t235 = t129 * t91 + t131 * t254;
t238 = (qJ(5) * t130 - t194) * qJD(1) + t235 + t248 * t200;
t220 = qJ(5) * t131;
t157 = qJ(6) * t129 - t220;
t145 = t157 * t132;
t237 = qJD(1) * t145 - qJD(4) * t157 + qJD(6) * t131 + t251;
t236 = qJ(5) * t199 - t205 * t220 + t251;
t234 = t129 * t94 + t78 * t199;
t90 = t115 * t213;
t233 = t129 * t78 + t90;
t232 = qJ(5) * t48;
t231 = qJ(5) * t85;
t229 = t112 * t85;
t228 = t112 * t87;
t227 = t129 * t54;
t226 = t131 * t54;
t223 = t57 * t129;
t222 = t57 * t131;
t216 = t129 * t130;
t221 = pkin(4) * t216 + t130 * t115;
t219 = t112 * t131;
t218 = t115 * t129;
t217 = t115 * t131;
t214 = t130 * t131;
t133 = qJD(3) ^ 2;
t212 = t133 * t130;
t211 = t133 * t132;
t210 = -qJD(5) - t23;
t16 = t24 - t247;
t208 = -qJD(6) - t16;
t207 = -t132 ^ 2 + t124;
t98 = qJD(1) * t116;
t195 = -t47 * t216 + t258 * t87;
t192 = pkin(8) * t112 * t129;
t191 = pkin(8) * t219;
t190 = pkin(8) * t202;
t189 = pkin(8) * t198;
t89 = t115 * t215;
t183 = t115 * t200;
t178 = -qJ(5) * t129 - pkin(3);
t177 = -pkin(4) - t218;
t176 = -t117 + t241;
t174 = t131 * t78 - t89;
t172 = -t129 * t77 - t131 * t56 - t58 * t199 + t55 * t200;
t169 = -qJ(5) + t217;
t167 = t258 * pkin(4) + qJ(5) * t182 + t115 * t201;
t59 = t85 * t184;
t164 = -qJD(4) * t90 + t131 * t94 - t78 * t200;
t33 = qJ(5) * t132 - t233;
t163 = -t102 - t172;
t11 = qJD(6) - t18 - t247;
t161 = t11 * t131 + t129 * t9;
t160 = -t11 * t129 + t131 * t9;
t17 = pkin(4) * t112 - t210;
t159 = t129 * t18 + t131 * t17;
t158 = t129 * t17 - t131 * t18;
t155 = 0.2e1 * qJD(3) * t98;
t152 = -pkin(5) * t48 - t172;
t151 = -pkin(5) * t47 + t171;
t3 = qJD(6) * t85 + t240 * t48 + t144;
t150 = t129 * t3 + t14 * t199;
t149 = -t131 * t3 + t14 * t200;
t148 = -t112 * t24 - t171;
t147 = t182 * t85 - t48 * t214 - t59;
t141 = -t14 * t85 + t152;
t26 = -t47 - t229;
t4 = -t111 - t163;
t139 = qJD(4) * t159 + t129 * t5 - t131 * t4;
t138 = t151 - t259;
t136 = -t228 - t48;
t134 = qJD(1) ^ 2;
t123 = t132 * pkin(4);
t109 = 0.2e1 * t111;
t100 = t248 * t129;
t96 = -pkin(4) * t131 + t178;
t76 = -t240 * t131 + t178;
t49 = -qJ(5) * t214 + t221;
t38 = pkin(4) * t87 + t231;
t35 = t130 * t157 + t221;
t34 = t123 - t174;
t31 = -pkin(5) * t216 - t33;
t30 = -pkin(4) * t206 - t175;
t29 = -qJ(5) * t206 - t235;
t27 = t240 * t87 + t231;
t25 = qJ(6) * t132 + t123 + t89 + (pkin(5) * t130 - t78) * t131;
t20 = (-qJ(5) * t201 - qJD(5) * t130) * t131 + t167;
t13 = t177 * t202 - t164;
t12 = (qJD(5) + t183) * t132 + t169 * t202 - t234;
t10 = qJD(3) * t145 + (qJD(6) * t129 + (qJ(6) * qJD(4) - qJD(5)) * t131) * t130 + t167;
t7 = -qJD(5) * t132 + (-pkin(5) * t214 - t89) * qJD(4) + (-t130 * t169 - t194) * qJD(3) + t234;
t6 = -pkin(5) * t182 + qJD(6) * t132 + (t193 + (-qJ(6) + t177) * t130) * qJD(3) - t164;
t2 = t152 + t253;
t1 = qJD(6) * t112 + t138;
t15 = [0, 0, 0, 0, 0.2e1 * t132 * t117, -0.2e1 * t207 * t197, t211, -t212, 0, -t115 * t211 + t130 * t155, t115 * t212 + t132 * t155, t87 * t184 + (-t131 * t47 - t87 * t200) * t130, t147 - t195, t146 - t252, t137 + t173 (-t112 - t205) * t202, -t164 * t112 + ((t115 * t85 + t227) * qJD(3) + t171) * t132 + (t54 * t199 + t115 * t48 + t223 + (qJD(1) * t174 - t112 * t218 - t23) * qJD(3)) * t130, t234 * t112 + (-t112 * t183 + (t115 * t87 + t226) * qJD(3) - t172) * t132 + (-t54 * t200 - t115 * t47 + t222 + (-t233 * qJD(1) - t112 * t217 - t24) * qJD(3)) * t130, t12 * t85 + t13 * t87 + t33 * t48 - t34 * t47 + t159 * t201 + (-qJD(4) * t158 + t129 * t4 + t131 * t5) * t130, -t112 * t13 - t20 * t85 - t48 * t49 + (-t22 * t203 - t5) * t132 + (-t22 * t199 - t245 + (qJD(1) * t34 + t17) * qJD(3)) * t130, t112 * t12 - t20 * t87 + t47 * t49 + (-t22 * t198 + t4) * t132 + (t22 * t200 - t244 + (-qJD(1) * t33 - t18) * qJD(3)) * t130, t12 * t18 + t13 * t17 + t20 * t22 + t33 * t4 + t34 * t5 + t49 * t8, -t25 * t47 - t31 * t48 + t6 * t87 - t7 * t85 + t160 * t201 + (-qJD(4) * t161 + t1 * t131 - t129 * t2) * t130, -t10 * t87 - t112 * t7 + t35 * t47 + (-t14 * t198 - t2) * t132 + ((qJD(1) * t31 + t11) * qJD(3) + t149) * t130, t10 * t85 + t112 * t6 + t35 * t48 + (t14 * t203 + t1) * t132 + ((-qJD(1) * t25 - t9) * qJD(3) + t150) * t130, t1 * t25 + t10 * t14 + t11 * t7 + t2 * t31 + t3 * t35 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, -t211, 0, 0, 0, 0, 0, t135, -t260, -t59 + (-t131 * t48 + t85 * t200) * t130 + t195, -t135, t260 (qJD(3) * t158 - t8) * t132 + (qJD(3) * t22 + t139) * t130, t147 + t195, t260, t135 (qJD(3) * t161 - t3) * t132 + (qJD(3) * t14 + qJD(4) * t160 + t1 * t129 + t131 * t2) * t130; 0, 0, 0, 0, -t130 * t134 * t132, t207 * t134, 0, 0, 0, qJD(3) * t65 - t98 * t206 - t57, -t98 * t205, -t129 * t47 - t87 * t219 (-t47 + t229) * t131 + (-t48 + t228) * t129, -t112 * t199 + (t112 * t213 + (-t87 + t203) * t130) * qJD(1), t112 * t200 + (-t112 * t215 + (t85 + t198) * t130) * qJD(1), t112 * t206, -pkin(3) * t48 - t222 + t175 * t112 - t65 * t85 + (t191 + t227) * qJD(4) + (t23 * t130 + (-t132 * t54 - t190) * t129) * qJD(1), pkin(3) * t47 + t223 - t235 * t112 - t65 * t87 + (-t192 + t226) * qJD(4) + (-t54 * t213 + (t24 - t189) * t130) * qJD(1), -t29 * t85 - t30 * t87 + (-t4 - t112 * t17 + (qJD(4) * t87 - t48) * pkin(8)) * t131 + ((qJD(4) * t85 - t47) * pkin(8) + t256) * t129, t112 * t30 + t244 - t48 * t96 + t236 * t85 + (-t129 * t22 - t191) * qJD(4) + (-t130 * t17 + (t132 * t22 + t190) * t129) * qJD(1), -t112 * t29 - t245 + t47 * t96 + t236 * t87 + (-t131 * t22 + t192) * qJD(4) + (t22 * t213 + (t18 + t189) * t130) * qJD(1), t139 * pkin(8) - t17 * t30 - t18 * t29 - t236 * t22 + t8 * t96, -t100 * t47 - t101 * t48 + t239 * t87 + t238 * t85 + (t2 - t246) * t131 + (t11 * t112 + t1) * t129, t47 * t76 + t237 * t87 + t238 * t112 + (t14 * t213 + (qJD(3) * t101 - t11) * t130) * qJD(1) - t150, t48 * t76 - t237 * t85 + t239 * t112 + (-t14 * t215 + (-qJD(3) * t100 + t9) * t130) * qJD(1) + t149, t1 * t100 + t101 * t2 - t238 * t11 - t237 * t14 + t239 * t9 + t3 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t82 - t249, t26, t136, t117, -t54 * t87 + t148, t112 * t23 + t54 * t85 + t172, pkin(4) * t47 - t232 + (-t18 - t24) * t87 + (t17 + t210) * t85, t38 * t85 - t148 - 0.2e1 * t166 + t242, t210 * t112 - t22 * t85 + t38 * t87 + t109 + t163, -pkin(4) * t5 - qJ(5) * t4 - t17 * t24 + t210 * t18 - t22 * t38, -t232 + t240 * t47 + (t11 + t208) * t87 + (t9 - t209) * t85, -t112 * t153 + t27 * t87 - 0.2e1 * t102 + t109 + t141, -t243 - t27 * t85 + (-0.2e1 * qJD(6) - t16) * t112 + 0.2e1 * t259 - t151, qJ(5) * t2 - t1 * t240 + t209 * t11 - t14 * t27 + t208 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t176, t255, t242 + t256, t26, t255, t176, t243 + (qJD(6) + t11) * t112 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t117 + t241, -t110 - t249, t141 - t246 + t253;];
tauc_reg  = t15;
