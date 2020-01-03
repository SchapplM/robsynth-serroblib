% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:48
% EndTime: 2019-12-31 19:10:55
% DurationCPUTime: 3.14s
% Computational Cost: add. (3455->349), mult. (8284->462), div. (0->0), fcn. (6564->14), ass. (0->188)
t151 = cos(pkin(9));
t158 = cos(qJ(3));
t211 = t158 * t151;
t132 = qJD(1) * t211;
t150 = sin(pkin(9));
t154 = sin(qJ(3));
t219 = t150 * t154;
t197 = qJD(1) * t219;
t108 = t132 - t197;
t248 = qJD(4) + qJD(5);
t262 = t108 - t248;
t100 = qJD(4) - t108;
t148 = pkin(9) + qJ(3);
t139 = sin(t148);
t140 = cos(t148);
t155 = sin(qJ(1));
t159 = cos(qJ(1));
t184 = g(1) * t159 + g(2) * t155;
t166 = -g(3) * t140 + t139 * t184;
t242 = pkin(6) + qJ(2);
t124 = t242 * t150;
t114 = qJD(1) * t124;
t125 = t242 * t151;
t115 = qJD(1) * t125;
t76 = -t154 * t114 + t158 * t115;
t259 = qJD(3) * t76;
t204 = qJD(1) * qJD(2);
t246 = qJDD(1) * t242 + t204;
t95 = t246 * t150;
t96 = t246 * t151;
t177 = -t154 * t96 - t158 * t95 - t259;
t27 = -qJDD(3) * pkin(3) - t177;
t261 = -qJD(4) * pkin(7) * t100 + t166 - t27;
t152 = sin(qJ(5));
t153 = sin(qJ(4));
t156 = cos(qJ(5));
t157 = cos(qJ(4));
t117 = t152 * t157 + t153 * t156;
t238 = t262 * t117;
t113 = t150 * t158 + t151 * t154;
t109 = t113 * qJD(1);
t205 = t157 * qJD(3);
t85 = t109 * t153 - t205;
t87 = qJD(3) * t153 + t109 * t157;
t179 = t152 * t85 - t156 * t87;
t37 = t152 * t87 + t156 * t85;
t260 = t179 * t37;
t209 = qJD(4) * t153;
t218 = t153 * t108;
t258 = t209 - t218;
t257 = t179 ^ 2 - t37 ^ 2;
t206 = qJD(5) * t156;
t207 = qJD(5) * t152;
t202 = t151 * qJDD(1);
t203 = t150 * qJDD(1);
t200 = qJD(3) * t132 + t154 * t202 + t158 * t203;
t73 = -qJD(3) * t197 + t200;
t33 = qJD(4) * t205 + qJDD(3) * t153 - t109 * t209 + t157 * t73;
t34 = qJD(4) * t87 - qJDD(3) * t157 + t153 * t73;
t7 = -t152 * t34 + t156 * t33 - t206 * t85 - t207 * t87;
t97 = qJD(5) + t100;
t256 = t37 * t97 + t7;
t134 = -pkin(2) * t151 - pkin(1);
t121 = qJD(1) * t134 + qJD(2);
t47 = -pkin(3) * t108 - pkin(7) * t109 + t121;
t69 = qJD(3) * pkin(7) + t76;
t21 = t153 * t47 + t157 * t69;
t16 = -pkin(8) * t85 + t21;
t14 = t16 * t207;
t149 = qJ(4) + qJ(5);
t145 = cos(t149);
t244 = g(3) * t139;
t251 = -t114 * t158 - t115 * t154;
t68 = -qJD(3) * pkin(3) - t251;
t35 = pkin(4) * t85 + t68;
t221 = t145 * t155;
t144 = sin(t149);
t222 = t144 * t159;
t90 = -t140 * t221 + t222;
t220 = t145 * t159;
t223 = t144 * t155;
t92 = t140 * t220 + t223;
t255 = g(1) * t92 - g(2) * t90 + t145 * t244 + t35 * t37 + t14;
t120 = qJDD(1) * t134 + qJDD(2);
t111 = t113 * qJD(3);
t180 = t154 * t203 - t158 * t202;
t74 = qJD(1) * t111 + t180;
t29 = pkin(3) * t74 - pkin(7) * t73 + t120;
t25 = t157 * t29;
t178 = -t154 * t95 + t158 * t96;
t26 = qJDD(3) * pkin(7) + qJD(3) * t251 + t178;
t67 = qJDD(4) + t74;
t2 = pkin(4) * t67 - pkin(8) * t33 - qJD(4) * t21 - t153 * t26 + t25;
t208 = qJD(4) * t157;
t172 = t153 * t29 + t157 * t26 + t208 * t47 - t209 * t69;
t3 = -pkin(8) * t34 + t172;
t199 = -t152 * t3 + t156 * t2;
t20 = -t153 * t69 + t157 * t47;
t15 = -pkin(8) * t87 + t20;
t11 = pkin(4) * t100 + t15;
t228 = t156 * t16;
t5 = t11 * t152 + t228;
t89 = t140 * t223 + t220;
t91 = -t140 * t222 + t221;
t254 = -g(1) * t91 + g(2) * t89 - qJD(5) * t5 + t144 * t244 + t35 * t179 + t199;
t164 = qJD(5) * t179 - t152 * t33 - t156 * t34;
t253 = -t179 * t97 + t164;
t59 = t117 * t113;
t112 = -t211 + t219;
t110 = t112 * qJD(3);
t213 = t157 * t110;
t250 = -t113 * t209 - t213;
t249 = qJ(2) * qJDD(1);
t116 = t152 * t153 - t156 * t157;
t239 = t262 * t116;
t64 = qJDD(5) + t67;
t247 = -t117 * t64 - t239 * t97;
t245 = pkin(7) + pkin(8);
t70 = pkin(3) * t109 - pkin(7) * t108;
t241 = t153 * t70 + t157 * t251;
t72 = pkin(3) * t112 - pkin(7) * t113 + t134;
t84 = -t124 * t154 + t125 * t158;
t77 = t157 * t84;
t240 = t153 * t72 + t77;
t237 = t100 * t85;
t236 = t100 * t87;
t235 = t109 * t37;
t234 = t109 * t179;
t233 = t109 * t85;
t232 = t109 * t87;
t230 = t153 * t33;
t229 = t153 * t67;
t227 = qJDD(1) * pkin(1);
t226 = t113 * t153;
t225 = t113 * t157;
t217 = t153 * t110;
t216 = t153 * t155;
t215 = t153 * t159;
t214 = t155 * t157;
t212 = t157 * t159;
t210 = t150 ^ 2 + t151 ^ 2;
t198 = qJD(4) * t245;
t195 = t113 * t208;
t194 = qJD(5) * t11 + t3;
t192 = -qJD(4) * t47 - t26;
t190 = t210 * qJD(1) ^ 2;
t189 = t100 * t157;
t188 = -t116 * t64 + t238 * t97;
t187 = 0.2e1 * t210;
t186 = pkin(4) * t258 - t76;
t185 = -t208 * t69 + t25;
t183 = g(1) * t155 - g(2) * t159;
t126 = t245 * t153;
t182 = -pkin(8) * t218 + qJD(5) * t126 + t153 * t198 + t241;
t127 = t245 * t157;
t57 = t157 * t70;
t181 = pkin(4) * t109 + qJD(5) * t127 - t153 * t251 + t57 + (-pkin(8) * t108 + t198) * t157;
t175 = -t124 * t158 - t125 * t154;
t174 = -t100 * t258 + t157 * t67;
t48 = -qJD(2) * t112 + qJD(3) * t175;
t71 = pkin(3) * t111 + pkin(7) * t110;
t171 = t153 * t71 + t157 * t48 + t208 * t72 - t209 * t84;
t170 = -pkin(7) * t67 + t100 * t68;
t169 = t195 - t217;
t168 = -t183 - t227;
t138 = qJDD(2) - t227;
t167 = -t138 - t168;
t163 = t187 * t204 - t184;
t49 = qJD(2) * t113 + qJD(3) * t84;
t137 = -pkin(4) * t157 - pkin(3);
t104 = t140 * t212 + t216;
t103 = -t140 * t215 + t214;
t102 = -t140 * t214 + t215;
t101 = t140 * t216 + t212;
t62 = t157 * t72;
t60 = t116 * t113;
t58 = t157 * t71;
t50 = pkin(4) * t226 - t175;
t28 = pkin(4) * t169 + t49;
t23 = -pkin(8) * t226 + t240;
t18 = pkin(4) * t112 - pkin(8) * t225 - t153 * t84 + t62;
t13 = -t207 * t226 + (t225 * t248 - t217) * t156 + t250 * t152;
t12 = t116 * t110 - t248 * t59;
t10 = pkin(4) * t34 + t27;
t9 = -pkin(8) * t169 + t171;
t6 = pkin(8) * t213 + pkin(4) * t111 - t153 * t48 + t58 + (-t77 + (pkin(8) * t113 - t72) * t153) * qJD(4);
t4 = t11 * t156 - t152 * t16;
t1 = [qJDD(1), t183, t184, t167 * t151, -t167 * t150, t187 * t249 + t163, (-t138 + t183) * pkin(1) + (t210 * t249 + t163) * qJ(2), -t109 * t110 + t113 * t73, -t108 * t110 - t109 * t111 - t112 * t73 - t113 * t74, -qJD(3) * t110 + qJDD(3) * t113, -qJD(3) * t111 - qJDD(3) * t112, 0, -qJD(3) * t49 + qJDD(3) * t175 + t111 * t121 + t112 * t120 + t134 * t74 + t140 * t183, -qJD(3) * t48 - qJDD(3) * t84 - t110 * t121 + t113 * t120 + t134 * t73 - t139 * t183, -t87 * t213 + (t157 * t33 - t209 * t87) * t113, -(-t153 * t87 - t157 * t85) * t110 + (-t230 - t157 * t34 + (t153 * t85 - t157 * t87) * qJD(4)) * t113, t100 * t250 + t111 * t87 + t112 * t33 + t225 * t67, -t100 * t169 - t111 * t85 - t112 * t34 - t226 * t67, t100 * t111 + t112 * t67, (-t84 * t208 + t58) * t100 + t62 * t67 + t185 * t112 + t20 * t111 + t49 * t85 - t175 * t34 + t68 * t195 - g(1) * t102 - g(2) * t104 + ((-qJD(4) * t72 - t48) * t100 - t84 * t67 + t192 * t112 + t27 * t113 - t68 * t110) * t153, -t171 * t100 - t240 * t67 - t172 * t112 - t21 * t111 + t49 * t87 - t175 * t33 - t68 * t213 - g(1) * t101 - g(2) * t103 + (t27 * t157 - t209 * t68) * t113, -t12 * t179 - t60 * t7, -t12 * t37 + t13 * t179 - t164 * t60 - t59 * t7, -t111 * t179 + t112 * t7 + t12 * t97 - t60 * t64, -t111 * t37 + t112 * t164 - t13 * t97 - t59 * t64, t111 * t97 + t112 * t64, (-t152 * t9 + t156 * t6) * t97 + (-t152 * t23 + t156 * t18) * t64 + t199 * t112 + t4 * t111 + t28 * t37 - t50 * t164 + t10 * t59 + t35 * t13 - g(1) * t90 - g(2) * t92 + ((-t152 * t18 - t156 * t23) * t97 - t5 * t112) * qJD(5), -g(1) * t89 - g(2) * t91 - t10 * t60 - t5 * t111 + t14 * t112 + t35 * t12 - t28 * t179 + t50 * t7 + (-(-qJD(5) * t23 + t6) * t97 - t18 * t64 - t2 * t112) * t152 + (-(qJD(5) * t18 + t9) * t97 - t23 * t64 - t194 * t112) * t156; 0, 0, 0, -t202, t203, -t190, -qJ(2) * t190 + qJDD(2) + t168, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t109 + t180, (t108 - t197) * qJD(3) + t200, 0, 0, 0, 0, 0, t174 - t233, -t100 ^ 2 * t157 - t229 - t232, 0, 0, 0, 0, 0, t188 - t235, t234 + t247; 0, 0, 0, 0, 0, 0, 0, -t109 * t108, -t108 ^ 2 + t109 ^ 2, (-t108 - t197) * qJD(3) + t200, -t180, qJDD(3), -t109 * t121 + t166 + t177 + t259, -t108 * t121 + t140 * t184 - t178 + t244, t189 * t87 + t230, (t33 - t237) * t157 + (-t34 - t236) * t153, t100 * t189 + t229 - t232, t174 + t233, -t100 * t109, -pkin(3) * t34 - t57 * t100 - t20 * t109 - t76 * t85 + (t100 * t251 + t170) * t153 + t261 * t157, -pkin(3) * t33 + t241 * t100 + t21 * t109 - t153 * t261 + t170 * t157 - t76 * t87, t117 * t7 - t179 * t239, -t116 * t7 + t117 * t164 - t179 * t238 - t239 * t37, t234 - t247, t188 + t235, -t97 * t109, (-t126 * t156 - t127 * t152) * t64 - t137 * t164 + t10 * t116 - t4 * t109 + (t152 * t182 - t156 * t181) * t97 + t186 * t37 - t238 * t35 + t166 * t145, -(-t126 * t152 + t127 * t156) * t64 + t137 * t7 + t10 * t117 + t5 * t109 + (t152 * t181 + t156 * t182) * t97 - t186 * t179 + t239 * t35 - t166 * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t85, -t85 ^ 2 + t87 ^ 2, t33 + t237, t236 - t34, t67, -g(1) * t103 + g(2) * t101 + t100 * t21 - t68 * t87 + (t192 + t244) * t153 + t185, g(1) * t104 - g(2) * t102 + t100 * t20 + t157 * t244 + t68 * t85 - t172, -t260, t257, t256, t253, t64, -(-t15 * t152 - t228) * t97 + (t156 * t64 - t207 * t97 - t37 * t87) * pkin(4) + t254, (-t16 * t97 - t2) * t152 + (t15 * t97 - t194) * t156 + (-t152 * t64 + t179 * t87 - t206 * t97) * pkin(4) + t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t260, t257, t256, t253, t64, t5 * t97 + t254, -t152 * t2 - t156 * t194 + t4 * t97 + t255;];
tau_reg = t1;
