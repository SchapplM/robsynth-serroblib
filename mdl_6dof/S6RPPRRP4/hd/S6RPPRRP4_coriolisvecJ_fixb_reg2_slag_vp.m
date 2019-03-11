% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:19
% EndTime: 2019-03-09 02:06:28
% DurationCPUTime: 3.08s
% Computational Cost: add. (3818->358), mult. (7408->471), div. (0->0), fcn. (4237->6), ass. (0->197)
t113 = sin(qJ(5));
t115 = cos(qJ(5));
t116 = cos(qJ(4));
t199 = t116 * qJD(1);
t100 = qJD(5) + t199;
t114 = sin(qJ(4));
t202 = t113 * qJD(4);
t209 = qJD(1) * t115;
t82 = t114 * t209 - t202;
t226 = t82 * t100;
t200 = t115 * qJD(4);
t210 = qJD(1) * t114;
t81 = t113 * t210 + t200;
t227 = t81 * t100;
t196 = qJD(1) * qJD(4);
t172 = t115 * t196;
t224 = qJD(5) * t81;
t57 = -t116 * t172 + t224;
t204 = qJD(5) * t115;
t174 = t114 * t204;
t179 = t116 * t202;
t132 = t174 + t179;
t58 = t132 * qJD(1) - qJD(5) * t202;
t264 = (t58 + t226) * t113 + (t57 + t227) * t115;
t110 = sin(pkin(9));
t111 = cos(pkin(9));
t178 = t114 * t200;
t219 = t115 * t116;
t220 = t113 * t116;
t78 = t110 * t220 + t115 * t111;
t243 = -qJD(5) * t78 - t110 * t178 - (t110 * t113 + t111 * t219) * qJD(1);
t201 = t114 * qJD(3);
t211 = qJD(1) * t111;
t117 = -pkin(1) - pkin(2);
t96 = t117 * qJD(1) + qJD(2);
t73 = qJ(2) * t211 + t110 * t96;
t61 = -qJD(1) * pkin(7) + t73;
t46 = t116 * t61 + t201;
t41 = qJD(4) * pkin(8) + t46;
t157 = t116 * pkin(4) + t114 * pkin(8);
t105 = t110 * qJ(2);
t72 = -qJD(1) * t105 + t111 * t96;
t60 = qJD(1) * pkin(3) - t72;
t42 = t157 * qJD(1) + t60;
t14 = t113 * t42 + t115 * t41;
t205 = qJD(5) * t113;
t197 = qJD(1) * qJD(2);
t171 = t111 * t197;
t258 = t116 * qJD(3) - t114 * t61;
t33 = qJD(4) * t258 + t116 * t171;
t156 = -pkin(4) * t114 + pkin(8) * t116;
t77 = t110 * qJD(2) + qJD(4) * t156;
t62 = t77 * qJD(1);
t168 = t113 * t33 - t115 * t62 + t41 * t204 + t42 * t205;
t137 = t14 * t100 - t168;
t206 = qJD(4) * t116;
t180 = t110 * t206;
t183 = t110 * t219;
t221 = t113 * t111;
t79 = t183 - t221;
t261 = (t110 * t57 + (qJD(4) * t79 + t111 * t82) * qJD(1)) * t114 - t243 * t100 - t82 * t180;
t186 = t100 * t220;
t260 = (t114 * (-t81 + t200) + t186) * qJD(1) + t100 * t205;
t12 = t100 * qJ(6) + t14;
t170 = t114 * t196;
t99 = pkin(5) * t170;
t2 = t99 + t168;
t259 = -t12 * t100 + t2;
t257 = t58 - t226;
t108 = t114 ^ 2;
t175 = t114 * t205;
t139 = t100 * t175 + t108 * t172;
t185 = t100 * t219;
t230 = t57 * t116;
t239 = t114 * t82;
t253 = qJD(4) * (t185 + t239) - t139 + t230;
t34 = t46 * qJD(4) + t114 * t171;
t233 = t34 * t114;
t123 = -(t114 * t46 + t116 * t258) * qJD(4) + t33 * t116 + t233;
t208 = qJD(2) * t111;
t181 = t116 * t208;
t203 = qJD(5) * t116;
t207 = qJD(4) * t114;
t164 = t111 * t117 - t105;
t84 = pkin(3) - t164;
t63 = t157 + t84;
t215 = t111 * qJ(2) + t110 * t117;
t85 = -pkin(7) + t215;
t9 = -t113 * (qJD(5) * t63 - t207 * t85 + t181) - t115 * (t203 * t85 - t77);
t252 = t82 ^ 2;
t251 = pkin(8) * t82;
t250 = pkin(8) * t113;
t40 = -qJD(4) * pkin(4) - t258;
t18 = -t81 * pkin(5) + t82 * qJ(6) + t40;
t249 = t18 * t82;
t5 = -t58 * pkin(5) - t57 * qJ(6) + t82 * qJD(6) + t34;
t248 = t5 * t113;
t247 = t5 * t115;
t246 = t81 * t82;
t153 = pkin(5) * t113 - qJ(6) * t115;
t245 = t201 + (-t153 * qJD(1) + t61) * t116 - qJD(5) * t153 + t113 * qJD(6);
t86 = t156 * qJD(1);
t23 = t113 * t86 + t115 * t258;
t229 = t58 * t115;
t177 = t116 * t200;
t67 = t81 * t177;
t244 = t114 * t229 + t67;
t223 = t110 * t114;
t242 = qJD(5) * t183 + t110 * t209 - t111 * t205 - t199 * t221 - t202 * t223;
t31 = t113 * t63 + t85 * t219;
t240 = t113 * t40;
t238 = t115 * t40;
t237 = t115 * t63;
t234 = t34 * t113;
t232 = t34 * t115;
t231 = t57 * t113;
t228 = t58 * t116;
t225 = qJD(1) * t84;
t119 = qJD(1) ^ 2;
t222 = t111 * t119;
t118 = qJD(4) ^ 2;
t218 = t118 * t114;
t217 = t118 * t116;
t13 = -t113 * t41 + t115 * t42;
t216 = qJD(6) - t13;
t109 = t116 ^ 2;
t214 = t108 - t109;
t213 = t108 + t109;
t212 = t118 + t119;
t195 = t113 * t77 + t115 * t181 + t63 * t204;
t69 = t82 * t174;
t194 = -t114 * t231 + t82 * t179 + t69;
t193 = t100 * t250;
t192 = pkin(8) * t100 * t115;
t191 = pkin(8) * t200;
t190 = 0.2e1 * t197;
t189 = 0.2e1 * t196;
t188 = t81 ^ 2 - t252;
t187 = t40 * t206;
t184 = t114 * t119 * t116;
t182 = t111 * t210;
t173 = t100 * t204;
t169 = t81 * t207 - t228;
t167 = t57 - t224;
t166 = -t60 - t208;
t163 = t111 * t189;
t162 = t110 * t190;
t160 = qJ(6) * t170;
t159 = t116 * t170;
t155 = t167 * pkin(8);
t154 = -t115 * pkin(5) - t113 * qJ(6);
t10 = -t100 * pkin(5) + t216;
t152 = t10 * t115 - t113 * t12;
t151 = t10 * t113 + t115 * t12;
t150 = t72 * t110 - t73 * t111;
t149 = t113 * t14 + t115 * t13;
t148 = -t113 * t13 + t115 * t14;
t22 = -t113 * t258 + t115 * t86;
t145 = t114 * t258 - t116 * t46;
t142 = qJD(1) * (t60 - t208);
t138 = -t153 + t85;
t135 = -t113 * t62 - t115 * t33 - t42 * t204 + t205 * t41;
t134 = -t118 * t85 + t162;
t133 = -t175 * t81 + t244;
t131 = qJD(4) * (t166 - t225);
t129 = -t113 * t227 + t229;
t127 = -t242 * t82 + t243 * t81 + t78 * t57 + t79 * t58;
t93 = t113 * t108 * t196;
t126 = -t100 * t132 + t93;
t1 = t100 * qJD(6) - t135 - t160;
t125 = qJD(5) * t152 + t1 * t115 + t2 * t113;
t124 = -qJD(5) * t149 + t113 * t168 - t115 * t135;
t122 = t81 * t174 + (t114 * t58 + t206 * t81) * t113;
t121 = t126 - t169;
t120 = -t58 * t223 - t81 * t180 - t242 * t100 + (qJD(4) * t78 + t111 * t81) * t210;
t92 = t170 * t250;
t89 = -pkin(4) + t154;
t88 = t100 * t210;
t64 = (-t100 - t199) * t207;
t52 = pkin(8) * t229;
t49 = -t82 * pkin(5) - t81 * qJ(6);
t43 = t138 * t114;
t30 = -t220 * t85 + t237;
t29 = t57 - t227;
t25 = -t237 + (t113 * t85 - pkin(5)) * t116;
t24 = t116 * qJ(6) + t31;
t21 = t173 + (t185 + (-t82 - t202) * t114) * qJD(1);
t20 = pkin(5) * t210 - t22;
t19 = -qJ(6) * t210 + t23;
t17 = -t115 * t226 + t231;
t16 = t138 * t206 + (qJD(5) * t154 + qJD(6) * t115 + t208) * t114;
t15 = t82 * t177 + (-t57 * t115 - t205 * t82) * t114;
t11 = t230 + (-t185 + t239) * qJD(4) + t139;
t8 = (-t113 * t203 - t178) * t85 + t195;
t7 = pkin(5) * t207 - t9;
t6 = (-t205 * t85 + qJD(6)) * t116 + (-t115 * t85 - qJ(6)) * t207 + t195;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, qJ(2) * t190, 0, 0, 0, 0, 0, 0, t162, 0.2e1 * t171, 0 ((-t110 * t164 + t111 * t215) * qJD(1) - t150) * qJD(2), 0.2e1 * t159, -t214 * t189, -t217, -0.2e1 * t159, t218, 0, t114 * t131 + t116 * t134, -t114 * t134 + t116 * t131, -t171 * t213 - t123, t123 * t85 + (-t145 * t111 + (t60 + t225) * t110) * qJD(2), t15, -t133 - t194, t11, t122, t114 * t173 + t228 - t93 + (-t114 * t81 + t186) * qJD(4), t64, t9 * t100 + (-t168 + (-t81 * t85 - t240) * qJD(4)) * t116 + (-t81 * t208 - t40 * t204 - t234 - t58 * t85 + (-qJD(1) * t30 - t13) * qJD(4)) * t114, -t8 * t100 + (t135 + (-t82 * t85 - t238) * qJD(4)) * t116 + (-t82 * t208 + t40 * t205 - t232 + t57 * t85 + (qJD(1) * t31 + t14) * qJD(4)) * t114, -t30 * t57 + t31 * t58 + t8 * t81 + t9 * t82 + t149 * t206 + (qJD(5) * t148 - t113 * t135 - t115 * t168) * t114, t85 * t187 + t13 * t9 + t14 * t8 - t135 * t31 - t168 * t30 + (t208 * t40 + t34 * t85) * t114, t15, t11, t67 + (-t205 * t81 + t229) * t114 + t194, t64, t126 + t169, t122, -t7 * t100 - t16 * t81 - t43 * t58 + (-t18 * t202 - t2) * t116 + (-t18 * t204 - t248 + (qJD(1) * t25 + t10) * qJD(4)) * t114, t24 * t58 + t25 * t57 + t6 * t81 - t7 * t82 - t152 * t206 + (qJD(5) * t151 + t1 * t113 - t115 * t2) * t114, t6 * t100 + t16 * t82 - t43 * t57 + (t18 * t200 + t1) * t116 + (-t18 * t205 + t247 + (-qJD(1) * t24 - t12) * qJD(4)) * t114, t1 * t24 + t10 * t7 + t12 * t6 + t18 * t16 + t2 * t25 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t119 * qJ(2), 0, 0, 0, 0, 0, 0, -t110 * t119, -t222, 0, t150 * qJD(1), 0, 0, 0, 0, 0, 0, -t110 * t116 * t212 + t114 * t163, t116 * t163 + t212 * t223, t213 * t222, t145 * t211 + (qJD(1) * t166 + t123) * t110, 0, 0, 0, 0, 0, 0, t120, t261, t127, -t40 * t182 - t135 * t79 + t168 * t78 + t243 * t14 - t242 * t13 + (t187 + t233) * t110, 0, 0, 0, 0, 0, 0, t120, t127, -t261, -t18 * t182 + t1 * t79 + t2 * t78 + t243 * t12 + (t114 * t5 + t18 * t206) * t110 + t242 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, -t217, 0, -qJD(4) * t145 + t33 * t114 - t34 * t116, 0, 0, 0, 0, 0, 0, t121, -t253, -t69 + (t114 * t167 - t206 * t82) * t113 + t244 (qJD(4) * t148 - t34) * t116 + (qJD(4) * t40 + t124) * t114, 0, 0, 0, 0, 0, 0, t121, t133 - t194, t253 (qJD(4) * t151 - t5) * t116 + (qJD(4) * t18 + t125) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t214 * t119, 0, t184, 0, 0, t114 * t142, t116 * t142, 0, 0, t17, t264, t21, t129, -t260, t88, pkin(4) * t58 - t22 * t100 - t232 + t46 * t81 + t92 + (-t192 + t240) * qJD(5) + (t114 * t13 + t220 * t40) * qJD(1), -pkin(4) * t57 + t23 * t100 + t234 + t46 * t82 + (t193 + t238) * qJD(5) + (t40 * t219 + (-t14 + t191) * t114) * qJD(1), -t22 * t82 - t23 * t81 + t52 + (-t13 * t199 - t135 + (-t13 - t251) * qJD(5)) * t115 + (t155 - t137) * t113, -t34 * pkin(4) + pkin(8) * t124 - t13 * t22 - t14 * t23 - t40 * t46, t17, t21, -t264, t88, t260, t129, t20 * t100 - t247 - t89 * t58 + t92 + t245 * t81 + (t113 * t18 - t192) * qJD(5) + (-t10 * t114 + t18 * t220) * qJD(1), -t19 * t81 + t20 * t82 + t52 + (t10 * t199 + t1 + (t10 - t251) * qJD(5)) * t115 + (t155 + t259) * t113, -t19 * t100 - t248 - t89 * t57 - t245 * t82 + (-t115 * t18 - t193) * qJD(5) + (-t18 * t219 + (t12 - t191) * t114) * qJD(1), pkin(8) * t125 - t10 * t20 - t12 * t19 - t18 * t245 + t5 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t188, t29, -t246, t257, -t170, t40 * t82 + t137, t13 * t100 - t40 * t81 + t135, 0, 0, t246, t29, t188, -t170, -t257, -t246, t49 * t81 + t137 + t249 - 0.2e1 * t99, -pkin(5) * t57 + t58 * qJ(6) + (-t12 + t14) * t82 + (-t10 + t216) * t81, -0.2e1 * t160 + t18 * t81 - t49 * t82 + (0.2e1 * qJD(6) - t13) * t100 - t135, -t2 * pkin(5) + t1 * qJ(6) - t10 * t14 + t12 * t216 - t18 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170 + t246, t29, -t100 ^ 2 - t252, -t249 + t259;];
tauc_reg  = t3;
