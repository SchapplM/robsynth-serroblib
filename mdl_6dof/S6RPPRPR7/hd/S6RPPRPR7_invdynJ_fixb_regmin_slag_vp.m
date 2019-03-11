% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:47
% EndTime: 2019-03-09 01:53:53
% DurationCPUTime: 3.17s
% Computational Cost: add. (4023->381), mult. (8238->492), div. (0->0), fcn. (6164->14), ass. (0->189)
t153 = pkin(9) + qJ(4);
t141 = sin(t153);
t143 = cos(t153);
t165 = sin(qJ(1));
t168 = cos(qJ(1));
t253 = g(1) * t165 - g(2) * t168;
t175 = -g(3) * t141 + t143 * t253;
t160 = cos(pkin(9));
t162 = -pkin(1) - qJ(3);
t125 = qJD(1) * t162 + qJD(2);
t204 = -pkin(7) * qJD(1) + t125;
t100 = t204 * t160;
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t215 = qJD(4) * t167;
t216 = qJD(4) * t164;
t158 = sin(pkin(9));
t247 = -qJD(1) * qJD(3) + qJDD(1) * t162;
t117 = qJDD(2) + t247;
t202 = -pkin(7) * qJDD(1) + t117;
t89 = t202 * t158;
t90 = t202 * t160;
t99 = t204 * t158;
t185 = -t100 * t216 - t164 * t89 + t167 * t90 - t215 * t99;
t25 = -qJDD(4) * pkin(4) + qJDD(5) - t185;
t258 = t25 + t175;
t115 = t158 * t167 + t160 * t164;
t104 = t115 * qJD(1);
t102 = qJD(6) + t104;
t163 = sin(qJ(6));
t166 = cos(qJ(6));
t225 = t167 * t160;
t207 = qJD(1) * t225;
t217 = qJD(1) * t158;
t208 = t164 * t217;
t106 = t207 - t208;
t157 = sin(pkin(10));
t159 = cos(pkin(10));
t83 = -qJD(4) * t159 + t106 * t157;
t85 = qJD(4) * t157 + t106 * t159;
t38 = t163 * t85 + t166 * t83;
t260 = t38 * t102;
t116 = t157 * t166 + t159 * t163;
t108 = t116 * qJD(6);
t234 = t104 * t116 + t108;
t188 = t163 * t83 - t166 * t85;
t259 = t102 * t188;
t109 = -t158 * t215 - t160 * t216;
t114 = t158 * t164 - t225;
t257 = t100 * t167 - t164 * t99;
t49 = -qJD(4) * pkin(4) + qJD(5) - t257;
t172 = t49 * t109 - t25 * t114 + t253;
t180 = t102 * t116;
t218 = t158 ^ 2 + t160 ^ 2;
t256 = t125 * t218;
t199 = g(1) * t168 + g(2) * t165;
t183 = t199 * t141;
t255 = t199 * t143;
t154 = qJDD(1) * qJ(2);
t155 = qJD(1) * qJD(2);
t252 = t154 + t155;
t123 = qJDD(3) + t252;
t254 = t123 - t199;
t251 = t104 * qJD(4);
t113 = t157 * t163 - t159 * t166;
t250 = t113 * qJD(6);
t249 = -qJD(6) + t102;
t235 = -t104 * t113 - t250;
t176 = -qJD(4) * t208 + qJDD(1) * t115;
t75 = qJD(4) * t207 + t176;
t73 = qJDD(6) + t75;
t239 = t116 * t73;
t248 = -t102 * t235 - t239;
t243 = g(3) * t143;
t173 = -t141 * t253 - t243;
t210 = t160 * qJDD(1);
t211 = t158 * qJDD(1);
t193 = -t164 * t211 + t167 * t210;
t74 = t193 - t251;
t63 = -qJDD(4) * t159 + t157 * t74;
t64 = qJDD(4) * t157 + t159 * t74;
t9 = -qJD(6) * t188 + t163 * t64 + t166 * t63;
t103 = t104 ^ 2;
t246 = 0.2e1 * t155;
t245 = pkin(8) * t159;
t144 = t158 * pkin(3);
t242 = -pkin(7) + t162;
t241 = pkin(8) + qJ(5);
t187 = t164 * t90 + t167 * t89;
t23 = qJDD(4) * qJ(5) + (qJD(5) + t257) * qJD(4) + t187;
t112 = pkin(3) * t211 + t123;
t24 = pkin(4) * t75 - qJ(5) * t74 - qJD(5) * t106 + t112;
t7 = t157 * t24 + t159 * t23;
t110 = -t158 * t216 + t160 * t215;
t43 = pkin(4) * t110 - qJ(5) * t109 + qJD(5) * t114 + qJD(2);
t118 = t242 * t158;
t119 = t242 * t160;
t78 = t118 * t164 - t119 * t167;
t50 = -qJD(3) * t115 - qJD(4) * t78;
t18 = t157 * t43 + t159 * t50;
t62 = t100 * t164 + t167 * t99;
t54 = qJD(4) * qJ(5) + t62;
t137 = qJD(1) * qJ(2) + qJD(3);
t120 = pkin(3) * t217 + t137;
t58 = pkin(4) * t104 - qJ(5) * t106 + t120;
t27 = t157 * t58 + t159 * t54;
t72 = pkin(4) * t106 + qJ(5) * t104;
t30 = t157 * t72 + t159 * t257;
t131 = qJ(2) + t144;
t71 = pkin(4) * t115 + qJ(5) * t114 + t131;
t79 = t118 * t167 + t119 * t164;
t35 = t157 * t71 + t159 * t79;
t240 = t106 * t38;
t48 = t113 * t73;
t238 = t157 * t75;
t237 = t159 * t75;
t236 = t188 * t106;
t233 = pkin(1) * qJDD(1);
t232 = t104 * t157;
t231 = t109 * t157;
t230 = t114 * t157;
t229 = t114 * t159;
t152 = pkin(10) + qJ(6);
t140 = sin(t152);
t228 = t165 * t140;
t142 = cos(t152);
t227 = t165 * t142;
t224 = t168 * t140;
t223 = t168 * t142;
t222 = -qJD(5) + t49;
t221 = qJD(4) * t109 - qJDD(4) * t114;
t220 = pkin(1) * t168 + qJ(2) * t165;
t214 = qJD(6) * t163;
t213 = qJD(6) * t166;
t6 = -t157 * t23 + t159 * t24;
t2 = pkin(5) * t75 - pkin(8) * t64 + t6;
t5 = -pkin(8) * t63 + t7;
t209 = -t163 * t5 + t166 * t2;
t206 = g(2) * t220;
t17 = -t157 * t50 + t159 * t43;
t26 = -t157 * t54 + t159 * t58;
t29 = -t157 * t257 + t159 * t72;
t34 = -t157 * t79 + t159 * t71;
t203 = t218 * t117;
t201 = qJDD(2) - t233;
t200 = -t102 * t234 - t48;
t197 = t7 * t157 + t6 * t159;
t196 = t163 * t2 + t166 * t5;
t195 = t26 * t110 + t6 * t115;
t194 = -t27 * t110 - t7 * t115;
t12 = pkin(5) * t104 - pkin(8) * t85 + t26;
t14 = -pkin(8) * t83 + t27;
t3 = t12 * t166 - t14 * t163;
t4 = t12 * t163 + t14 * t166;
t191 = -t157 * t26 + t159 * t27;
t22 = pkin(5) * t115 + pkin(8) * t229 + t34;
t28 = pkin(8) * t230 + t35;
t190 = -t163 * t28 + t166 * t22;
t189 = t163 * t22 + t166 * t28;
t184 = -qJD(4) * t110 - qJDD(4) * t115;
t8 = -t163 * t63 + t166 * t64 - t213 * t83 - t214 * t85;
t121 = t241 * t157;
t182 = pkin(8) * t232 - qJD(5) * t159 + qJD(6) * t121 + t30;
t122 = t241 * t159;
t181 = pkin(5) * t106 + qJD(5) * t157 + qJD(6) * t122 + t104 * t245 + t29;
t179 = t102 * t113;
t177 = pkin(4) * t141 - qJ(5) * t143 + t144;
t170 = t252 + t254;
t51 = -qJD(3) * t114 + qJD(4) * t79;
t169 = qJD(1) ^ 2;
t161 = -pkin(7) - qJ(3);
t146 = t168 * qJ(2);
t134 = -pkin(5) * t159 - pkin(4);
t94 = t141 * t223 - t228;
t93 = t141 * t224 + t227;
t92 = t141 * t227 + t224;
t91 = -t141 * t228 + t223;
t70 = t113 * t114;
t69 = t116 * t114;
t52 = -pkin(5) * t230 + t78;
t37 = -pkin(5) * t232 + t62;
t36 = pkin(5) * t231 + t51;
t33 = pkin(5) * t83 + t49;
t32 = t109 * t116 - t213 * t229 + t214 * t230;
t31 = t108 * t114 - t109 * t113;
t13 = -pkin(8) * t231 + t18;
t11 = pkin(5) * t110 - t109 * t245 + t17;
t10 = pkin(5) * t63 + t25;
t1 = [qJDD(1), t253, t199, qJDD(2) - t253 - 0.2e1 * t233, 0.2e1 * t154 + t246 - t199, -t201 * pkin(1) - g(1) * (-t165 * pkin(1) + t146) - t206 + (t154 + t246) * qJ(2), t170 * t158, t170 * t160, t253 + t218 * (-t117 - t247) t123 * qJ(2) + t137 * qJD(2) - g(1) * (t162 * t165 + t146) - g(2) * (qJ(3) * t168 + t220) + t162 * t203 - qJD(3) * t256, t106 * t109 - t114 * t74, -t104 * t109 - t106 * t110 + t114 * t75 - t115 * t74, t221, t184, 0, qJD(2) * t104 - t51 * qJD(4) - t78 * qJDD(4) + t120 * t110 + t112 * t115 + t131 * t75 - t183, qJD(2) * t106 - t50 * qJD(4) - t79 * qJDD(4) + t120 * t109 - t112 * t114 + t131 * t74 - t255, t17 * t104 + t157 * t172 - t159 * t183 + t34 * t75 + t51 * t83 + t78 * t63 + t195, -t18 * t104 + t157 * t183 + t159 * t172 - t35 * t75 + t51 * t85 + t78 * t64 + t194, -t17 * t85 - t18 * t83 - t34 * t64 - t35 * t63 + t255 + t197 * t114 + (-t157 * t27 - t159 * t26) * t109, t7 * t35 + t27 * t18 + t6 * t34 + t26 * t17 + t25 * t78 + t49 * t51 - g(1) * t146 - t206 + (-g(1) * t177 + g(2) * t161) * t168 + (-g(1) * (-pkin(1) + t161) - g(2) * t177) * t165, -t188 * t31 + t70 * t8, t188 * t32 - t31 * t38 + t69 * t8 - t70 * t9, t102 * t31 - t110 * t188 + t115 * t8 + t70 * t73, -t102 * t32 - t110 * t38 - t115 * t9 + t69 * t73, t102 * t110 + t115 * t73 (t166 * t11 - t163 * t13) * t102 + t190 * t73 + t209 * t115 + t3 * t110 + t36 * t38 + t52 * t9 - t10 * t69 + t33 * t32 - g(1) * t94 - g(2) * t92 + (-t102 * t189 - t115 * t4) * qJD(6) -(t163 * t11 + t166 * t13) * t102 - t189 * t73 - t196 * t115 - t4 * t110 - t36 * t188 + t52 * t8 + t10 * t70 + t33 * t31 + g(1) * t93 - g(2) * t91 + (-t102 * t190 - t115 * t3) * qJD(6); 0, 0, 0, qJDD(1), -t169, -qJ(2) * t169 + t201 - t253, -t169 * t158, -t169 * t160, -t218 * qJDD(1), -t137 * qJD(1) + t203 - t253, 0, 0, 0, 0, 0, -qJD(1) * t104 + t221, -qJD(1) * t106 + t184, -t115 * t238 - t109 * t83 + t114 * t63 + (-qJD(1) * t159 - t110 * t157) * t104, -t115 * t237 - t109 * t85 + t114 * t64 + (qJD(1) * t157 - t110 * t159) * t104 (qJD(1) * t85 - t110 * t83 - t115 * t63) * t159 + (qJD(1) * t83 + t110 * t85 + t115 * t64) * t157 (-qJD(1) * t26 - t194) * t159 + (-qJD(1) * t27 - t195) * t157 - t172, 0, 0, 0, 0, 0, -t109 * t38 + t114 * t9 - t110 * t180 + qJD(1) * t179 + (t102 * t250 - t239) * t115, t109 * t188 + t114 * t8 + t110 * t179 + qJD(1) * t180 + (qJD(6) * t180 + t48) * t115; 0, 0, 0, 0, 0, 0, t211, t210, -t218 * t169, qJD(1) * t256 + t254, 0, 0, 0, 0, 0 (t106 + t207) * qJD(4) + t176, t193 - 0.2e1 * t251, -t103 * t157 - t106 * t83 + t237, -t103 * t159 - t106 * t85 - t238, -t157 * t63 - t159 * t64 + (t157 * t85 - t159 * t83) * t104, t104 * t191 - t49 * t106 + t197 - t199, 0, 0, 0, 0, 0, t200 - t240, t236 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t104, t106 ^ 2 - t103, t193 (t106 - t207) * qJD(4) - t176, qJDD(4), t62 * qJD(4) - t120 * t106 - t175 + t185, t120 * t104 - t173 - t187, -qJ(5) * t238 - pkin(4) * t63 - t26 * t106 - t62 * t83 + (t157 * t222 - t29) * t104 - t258 * t159, -qJ(5) * t237 - pkin(4) * t64 + t27 * t106 - t62 * t85 + (t159 * t222 + t30) * t104 + t258 * t157, t29 * t85 + t30 * t83 + (-qJ(5) * t63 - qJD(5) * t83 - t104 * t26 + t7) * t159 + (qJ(5) * t64 + qJD(5) * t85 - t104 * t27 - t6) * t157 + t173, -t26 * t29 - t27 * t30 - t49 * t62 + t191 * qJD(5) - t258 * pkin(4) + (-t6 * t157 + t7 * t159 + t173) * qJ(5), t8 * t116 - t188 * t235, -t8 * t113 - t116 * t9 + t188 * t234 - t235 * t38, t236 - t248, t200 + t240, -t102 * t106 (-t121 * t166 - t122 * t163) * t73 + t134 * t9 + t10 * t113 - t3 * t106 - t37 * t38 + t234 * t33 + (t163 * t182 - t166 * t181) * t102 - t175 * t142 -(-t121 * t163 + t122 * t166) * t73 + t134 * t8 + t10 * t116 + t4 * t106 + t37 * t188 + t235 * t33 + (t163 * t181 + t166 * t182) * t102 + t175 * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * t85 + t63, -t104 * t83 + t64, -t83 ^ 2 - t85 ^ 2, t26 * t85 + t27 * t83 + t258, 0, 0, 0, 0, 0, t9 - t259, t8 - t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188 * t38, t188 ^ 2 - t38 ^ 2, t8 + t260, -t9 - t259, t73, -g(1) * t91 - g(2) * t93 + t140 * t243 + t188 * t33 + t249 * t4 + t209, g(1) * t92 - g(2) * t94 + t142 * t243 + t249 * t3 + t33 * t38 - t196;];
tau_reg  = t1;
