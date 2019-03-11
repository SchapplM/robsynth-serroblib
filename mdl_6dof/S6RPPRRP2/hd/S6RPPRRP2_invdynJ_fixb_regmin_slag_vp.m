% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:20
% EndTime: 2019-03-09 02:01:27
% DurationCPUTime: 2.79s
% Computational Cost: add. (4707->386), mult. (9967->460), div. (0->0), fcn. (7310->14), ass. (0->194)
t138 = pkin(10) + qJ(4);
t130 = cos(t138);
t243 = g(3) * t130;
t128 = sin(t138);
t139 = qJ(1) + pkin(9);
t129 = sin(t139);
t131 = cos(t139);
t190 = g(1) * t131 + g(2) * t129;
t256 = t190 * t128;
t160 = -t243 + t256;
t146 = sin(qJ(4));
t149 = cos(qJ(4));
t141 = sin(pkin(9));
t118 = t141 * pkin(1) + qJ(3);
t112 = t118 * qJD(1);
t142 = cos(pkin(10));
t127 = t142 * qJD(2);
t140 = sin(pkin(10));
t77 = t127 + (-pkin(7) * qJD(1) - t112) * t140;
t209 = qJD(1) * t142;
t91 = t140 * qJD(2) + t142 * t112;
t78 = pkin(7) * t209 + t91;
t38 = t146 * t77 + t149 * t78;
t266 = t38 * qJD(4);
t105 = qJD(1) * qJD(3) + t118 * qJDD(1);
t125 = t142 * qJDD(2);
t71 = t125 + (-pkin(7) * qJDD(1) - t105) * t140;
t201 = t142 * qJDD(1);
t84 = t140 * qJDD(2) + t142 * t105;
t72 = pkin(7) * t201 + t84;
t177 = -t146 * t72 + t149 * t71 - t266;
t12 = -qJDD(4) * pkin(4) - t177;
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t213 = t146 * t140;
t107 = -t149 * t142 + t213;
t99 = t107 * qJD(4);
t164 = qJD(1) * t99;
t108 = t149 * t140 + t146 * t142;
t165 = t108 * qJDD(1);
t153 = t165 - t164;
t205 = t148 * qJD(4);
t208 = qJD(5) * t145;
t254 = t108 * qJD(1);
t30 = -qJD(5) * t205 - t145 * qJDD(4) - t148 * t153 + t208 * t254;
t166 = t107 * qJD(1);
t206 = t145 * qJD(4);
t207 = qJD(5) * t148;
t31 = -t148 * qJDD(4) + t145 * t165 + (qJD(5) - t166) * t206 + t207 * t254;
t81 = t148 * t254 + t206;
t3 = t31 * pkin(5) + t30 * qJ(6) - t81 * qJD(6) + t12;
t267 = -t3 + t160;
t143 = cos(pkin(9));
t240 = t143 * pkin(1);
t122 = -pkin(2) - t240;
t203 = qJDD(1) * t122;
t109 = qJDD(3) + t203;
t196 = -g(1) * t129 + g(2) * t131;
t265 = -t109 - t196;
t97 = -qJD(1) * t213 + t149 * t209;
t93 = qJD(5) - t97;
t264 = qJD(4) * t254;
t180 = t146 * t71 + t149 * t72;
t258 = -t146 * t78 + t149 * t77;
t11 = qJDD(4) * pkin(8) + qJD(4) * t258 + t180;
t121 = t142 * pkin(3) + pkin(2);
t110 = -t121 - t240;
t255 = -t108 * pkin(8) + t110;
t100 = t108 * qJD(4);
t202 = t140 * qJDD(1);
t183 = t146 * t202 - t149 * t201;
t65 = qJD(1) * t100 + t183;
t24 = t65 * pkin(4) + pkin(8) * t164 + t255 * qJDD(1) + qJDD(3);
t35 = qJD(4) * pkin(8) + t38;
t95 = t110 * qJD(1) + qJD(3);
t45 = -t97 * pkin(4) - pkin(8) * t254 + t95;
t197 = t145 * t11 - t148 * t24 + t35 * t207 + t45 * t208;
t62 = qJDD(5) + t65;
t249 = t62 * pkin(5);
t2 = qJDD(6) + t197 - t249;
t14 = t145 * t45 + t148 * t35;
t8 = t93 * qJ(6) + t14;
t248 = t8 * t93;
t262 = -t2 + t248;
t79 = t145 * t254 - t205;
t261 = -t100 * t79 - t107 * t31;
t53 = t145 * t62;
t260 = -t93 * t207 - t53;
t198 = t93 * t208;
t54 = t148 * t62;
t259 = t198 - t54;
t257 = t196 * t128;
t34 = -qJD(4) * pkin(4) - t258;
t17 = t79 * pkin(5) - t81 * qJ(6) + t34;
t250 = pkin(8) * t62;
t253 = t17 * t93 - t250;
t244 = g(3) * t128;
t161 = -t190 * t130 - t244;
t252 = t81 ^ 2;
t251 = t93 ^ 2;
t241 = t14 * t93;
t147 = sin(qJ(1));
t239 = t147 * pkin(1);
t238 = t79 * t97;
t237 = t81 * t79;
t195 = t81 * t93;
t236 = t81 * t254;
t235 = t254 * t79;
t234 = pkin(7) + t118;
t218 = t108 * t148;
t222 = t148 * t99;
t233 = -t31 * t218 + t79 * t222;
t232 = -t145 * t31 - t79 * t207;
t63 = pkin(4) * t254 - t97 * pkin(8);
t231 = t145 * t63 + t148 * t258;
t230 = t81 * t100 - t107 * t30;
t52 = t107 * pkin(4) + t255;
t101 = t234 * t140;
t102 = t234 * t142;
t58 = -t146 * t101 + t149 * t102;
t229 = t145 * t52 + t148 * t58;
t184 = pkin(5) * t145 - qJ(6) * t148;
t228 = -t145 * qJD(6) + t93 * t184 - t38;
t227 = pkin(8) * qJD(5);
t226 = t145 * t79;
t225 = t145 * t93;
t224 = t145 * t99;
t223 = t148 * t81;
t220 = t30 * t145;
t219 = t62 * qJ(6);
t217 = t129 * t145;
t216 = t129 * t148;
t215 = t131 * t145;
t214 = t131 * t148;
t13 = -t145 * t35 + t148 * t45;
t212 = qJD(6) - t13;
t211 = (g(1) * t214 + g(2) * t216) * t128;
t210 = t140 ^ 2 + t142 ^ 2;
t200 = t93 * t227;
t61 = t93 * t222;
t199 = t81 * t224;
t86 = t130 * t217 + t214;
t88 = t130 * t215 - t216;
t193 = -g(1) * t86 + g(2) * t88;
t87 = t130 * t216 - t215;
t89 = t130 * t214 + t217;
t192 = g(1) * t87 - g(2) * t89;
t170 = t148 * t11 + t145 * t24 + t45 * t207 - t35 * t208;
t1 = t93 * qJD(6) + t170 + t219;
t7 = -t93 * pkin(5) + t212;
t191 = t93 * t7 + t1;
t150 = cos(qJ(1));
t188 = g(1) * t147 - g(2) * t150;
t187 = -t145 * t8 + t148 * t7;
t186 = t145 * t7 + t148 * t8;
t185 = t148 * pkin(5) + t145 * qJ(6);
t83 = -t140 * t105 + t125;
t182 = -t83 * t140 + t84 * t142;
t181 = (-t140 * t112 + t127) * t140 - t91 * t142;
t178 = -t148 * t97 * t93 - t260;
t176 = -t149 * t101 - t146 * t102;
t175 = t97 * t225 - t259;
t174 = t200 + t243;
t173 = pkin(4) + t185;
t172 = t130 * pkin(4) + t128 * pkin(8) + t121;
t171 = t93 * t34 - t250;
t41 = -t107 * qJD(3) + t176 * qJD(4);
t64 = t100 * pkin(4) + t99 * pkin(8);
t169 = t145 * t64 + t148 * t41 + t52 * t207 - t58 * t208;
t168 = -t108 * t198 + t62 * t218 - t61;
t163 = -t203 + t265;
t162 = g(1) * t88 + g(2) * t86 + t145 * t244 - t197;
t159 = t260 * t108 + t93 * t224;
t158 = t187 * qJD(5) + t1 * t148 + t2 * t145;
t157 = t17 * t81 + qJDD(6) - t162;
t155 = t159 - t261;
t154 = -g(1) * t89 - g(2) * t87 - t148 * t244 + t170;
t42 = t108 * qJD(3) + t58 * qJD(4);
t144 = -pkin(7) - qJ(3);
t135 = t150 * pkin(1);
t92 = t110 * qJDD(1) + qJDD(3);
t67 = -t100 * qJD(4) - t107 * qJDD(4);
t66 = -t99 * qJD(4) + t108 * qJDD(4);
t43 = t81 * pkin(5) + t79 * qJ(6);
t29 = t184 * t108 - t176;
t20 = -t107 * pkin(5) + t145 * t58 - t148 * t52;
t19 = t107 * qJ(6) + t229;
t18 = t79 * t93 - t30;
t16 = -pkin(5) * t254 + t145 * t258 - t148 * t63;
t15 = qJ(6) * t254 + t231;
t6 = -t184 * t99 + (t185 * qJD(5) - qJD(6) * t148) * t108 + t42;
t5 = -t100 * pkin(5) + t229 * qJD(5) + t145 * t41 - t148 * t64;
t4 = t100 * qJ(6) + t107 * qJD(6) + t169;
t9 = [qJDD(1), t188, g(1) * t150 + g(2) * t147 (t188 + (t141 ^ 2 + t143 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t163 * t142, -t163 * t140, t105 * t210 + t182 - t190, t109 * t122 - g(1) * (-t129 * pkin(2) + t131 * qJ(3) - t239) - g(2) * (t131 * pkin(2) + t129 * qJ(3) + t135) + t182 * t118 - t181 * qJD(3), t108 * t153 - t254 * t99, -t100 * t254 - t107 * t153 - t108 * t65 - t99 * t97, t66, t67, 0, -t42 * qJD(4) + qJDD(4) * t176 + t95 * t100 + t92 * t107 + t110 * t65 - t130 * t196, -t58 * qJDD(4) + t92 * t108 - t95 * t99 + t257 + t110 * t165 + (-t110 * t166 - t41) * qJD(4), -t81 * t222 + (-t30 * t148 - t81 * t208) * t108, t199 + (t220 + (-t223 + t226) * qJD(5)) * t108 + t233, t168 + t230, t159 + t261, t93 * t100 + t62 * t107, -t197 * t107 + t13 * t100 + t42 * t79 - t176 * t31 + ((-qJD(5) * t58 + t64) * t93 + t52 * t62 + t34 * qJD(5) * t108) * t148 + ((-qJD(5) * t52 - t41) * t93 - t58 * t62 + t12 * t108 - t34 * t99) * t145 + t192, -t169 * t93 - t229 * t62 - t170 * t107 - t14 * t100 + t42 * t81 + t176 * t30 - t34 * t222 + (t12 * t148 - t208 * t34) * t108 + t193, -t17 * t224 - t7 * t100 - t2 * t107 - t20 * t62 + t29 * t31 - t5 * t93 + t6 * t79 + (t3 * t145 + t17 * t207) * t108 + t192, -t19 * t31 - t20 * t30 - t4 * t79 + t5 * t81 - t187 * t99 - t257 + (-qJD(5) * t186 - t1 * t145 + t2 * t148) * t108, t17 * t222 + t1 * t107 + t8 * t100 + t19 * t62 + t29 * t30 + t4 * t93 - t6 * t81 + (-t3 * t148 + t17 * t208) * t108 - t193, t1 * t19 + t8 * t4 + t3 * t29 + t17 * t6 + t2 * t20 + t7 * t5 - g(1) * (-t87 * pkin(5) - t86 * qJ(6) - t239) - g(2) * (t89 * pkin(5) + t88 * qJ(6) + t135) + (g(1) * t144 - g(2) * t172) * t131 + (g(1) * t172 + g(2) * t144) * t129; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t84 * t140 + t83 * t142 - g(3), 0, 0, 0, 0, 0, t67, -t66, 0, 0, 0, 0, 0, t155, t259 * t108 + t230 + t61, t155, -t199 + (-t220 + (t223 + t226) * qJD(5)) * t108 + t233, t168 - t230, t17 * t100 + t3 * t107 + t108 * t158 - t186 * t99 - g(3); 0, 0, 0, 0, -t201, t202, -t210 * qJD(1) ^ 2, t181 * qJD(1) - t265, 0, 0, 0, 0, 0, t183 + 0.2e1 * t264, t165 + (t97 - t166) * qJD(4), 0, 0, 0, 0, 0, t175 - t235, -t148 * t251 - t236 - t53, -t225 * t93 - t235 + t54 (t30 + t238) * t148 + t145 * t195 + t232, t178 + t236, t191 * t145 + t262 * t148 - t17 * t254 + t196; 0, 0, 0, 0, 0, 0, 0, 0, -t254 * t97, t254 ^ 2 - t97 ^ 2, t165 + (-t97 - t166) * qJD(4), -t183, qJDD(4), -t254 * t95 + t160 + t177 + t266, -t95 * t97 - t161 - t180, t148 * t195 - t220 (-t30 + t238) * t148 - t81 * t225 + t232, t178 - t236, t175 + t235, -t93 * t254, -pkin(4) * t31 - t13 * t254 - t38 * t79 + (-t243 - t12 + (-t63 - t227) * t93) * t148 + (t258 * t93 + t171) * t145 + t211, pkin(4) * t30 + t231 * t93 + t14 * t254 - t38 * t81 + t171 * t148 + (t12 + t174 - t256) * t145, -t173 * t31 + t16 * t93 + t7 * t254 + t228 * t79 + (-t174 - t3) * t148 + t253 * t145 + t211, t15 * t79 - t16 * t81 + ((qJD(5) * t81 - t31) * pkin(8) + t191) * t148 + ((qJD(5) * t79 - t30) * pkin(8) - t262) * t145 + t161, -t173 * t30 - t15 * t93 - t8 * t254 - t228 * t81 - t253 * t148 + (-t200 + t267) * t145, -t8 * t15 - t7 * t16 + t228 * t17 + (t158 + t161) * pkin(8) + t267 * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, -t79 ^ 2 + t252, t18, t195 - t31, t62, -t34 * t81 + t162 + t241, t13 * t93 + t34 * t79 - t154, -t43 * t79 - t157 + t241 + 0.2e1 * t249, pkin(5) * t30 - t31 * qJ(6) + (-t14 + t8) * t81 + (t7 - t212) * t79, 0.2e1 * t219 - t17 * t79 + t43 * t81 + (0.2e1 * qJD(6) - t13) * t93 + t154, t1 * qJ(6) - t2 * pkin(5) - t17 * t43 - t7 * t14 - g(1) * (-t88 * pkin(5) + t89 * qJ(6)) - g(2) * (-t86 * pkin(5) + t87 * qJ(6)) + t212 * t8 + t184 * t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t183 + t237 - t264, t18, -t251 - t252, t157 - t248 - t249;];
tau_reg  = t9;
