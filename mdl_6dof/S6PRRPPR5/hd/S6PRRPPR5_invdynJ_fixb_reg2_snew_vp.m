% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPPR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:27:43
% EndTime: 2019-05-05 03:27:54
% DurationCPUTime: 4.20s
% Computational Cost: add. (13291->359), mult. (28268->499), div. (0->0), fcn. (18504->12), ass. (0->217)
t202 = sin(qJ(3));
t188 = t202 ^ 2;
t208 = qJD(2) ^ 2;
t184 = t188 * t208;
t207 = qJD(3) ^ 2;
t172 = -t184 - t207;
t205 = cos(qJ(3));
t249 = t205 * t208;
t234 = t202 * t249;
t168 = qJDD(3) - t234;
t250 = t205 * t168;
t125 = t202 * t172 + t250;
t239 = qJD(2) * qJD(3);
t181 = t205 * t239;
t183 = t202 * qJDD(2);
t156 = t183 + 0.2e1 * t181;
t195 = sin(pkin(6));
t198 = cos(pkin(6));
t203 = sin(qJ(2));
t206 = cos(qJ(2));
t291 = t198 * (t202 * t168 - t205 * t172) + (t203 * t125 + t206 * t156) * t195;
t290 = pkin(8) * t125;
t189 = t205 ^ 2;
t185 = t189 * t208;
t174 = -t185 - t207;
t167 = qJDD(3) + t234;
t253 = t202 * t167;
t124 = -t205 * t174 + t253;
t231 = t202 * t239;
t237 = t205 * qJDD(2);
t159 = -0.2e1 * t231 + t237;
t288 = (t203 * t124 - t206 * t159) * t195 - t198 * (t205 * t167 + t202 * t174);
t287 = -2 * qJD(5);
t286 = pkin(8) * t124;
t193 = sin(pkin(11));
t196 = cos(pkin(11));
t244 = qJD(2) * t205;
t146 = t193 * qJD(3) + t196 * t244;
t148 = t196 * qJD(3) - t193 * t244;
t119 = t148 * t146;
t157 = t183 + t181;
t276 = -t119 + t157;
t285 = t193 * t276;
t284 = t196 * t276;
t201 = sin(qJ(6));
t149 = qJDD(6) + t157;
t204 = cos(qJ(6));
t113 = t204 * t146 + t201 * t148;
t115 = -t201 * t146 + t204 * t148;
t91 = t115 * t113;
t279 = -t91 + t149;
t283 = t201 * t279;
t282 = t204 * t279;
t194 = sin(pkin(10));
t197 = cos(pkin(10));
t164 = t194 * g(1) - t197 * g(2);
t190 = -g(3) + qJDD(1);
t281 = t164 * t198 + t190 * t195;
t255 = t202 * qJ(4);
t221 = -t205 * pkin(3) - t255;
t165 = -t197 * g(1) - t194 * g(2);
t110 = t206 * t165 + t203 * t281;
t96 = -t208 * pkin(2) + qJDD(2) * pkin(8) + t110;
t228 = t208 * t221 + t96;
t275 = -t207 * pkin(3) + t228 * t205;
t133 = -t195 * t164 + t198 * t190;
t127 = t205 * t133;
t226 = -qJDD(3) * pkin(3) - t207 * qJ(4) + qJDD(4) - t127;
t60 = -qJDD(3) * qJ(5) + (t157 - t181) * pkin(4) + (-qJ(5) * t249 + t228) * t202 + t226;
t158 = -t231 + t237;
t241 = t202 * qJD(2);
t166 = pkin(4) * t241 - qJD(3) * qJ(5);
t224 = t203 * t165 - t206 * t281;
t95 = -qJDD(2) * pkin(2) - t208 * pkin(8) + t224;
t211 = -t158 * pkin(3) + t95 + (-t157 - t181) * qJ(4);
t229 = pkin(3) * qJD(3) - (2 * qJD(4));
t69 = -pkin(4) * t185 - t158 * qJ(5) + (-t166 + t229) * t241 + t211;
t225 = t148 * t287 - t193 * t69 + t196 * t60;
t33 = t146 * t287 + t193 * t60 + t196 * t69;
t274 = t250 + t202 * (t185 - t207);
t111 = t113 ^ 2;
t112 = t115 ^ 2;
t143 = t146 ^ 2;
t144 = t148 ^ 2;
t177 = qJD(6) + t241;
t175 = t177 ^ 2;
t273 = 2 * qJD(4);
t131 = t196 * qJDD(3) - t193 * t158;
t233 = t146 * t241;
t222 = -t131 - t233;
t209 = t276 * pkin(5) + t222 * pkin(9) + t225;
t130 = -t193 * qJDD(3) - t196 * t158;
t132 = pkin(5) * t241 - t148 * pkin(9);
t26 = -t143 * pkin(5) + t130 * pkin(9) - t132 * t241 + t33;
t11 = t201 * t26 - t204 * t209;
t12 = t201 * t209 + t204 * t26;
t8 = -t204 * t11 + t201 * t12;
t272 = t193 * t8;
t271 = t196 * t8;
t270 = -pkin(3) - qJ(5);
t238 = qJDD(3) * qJ(4);
t254 = t202 * t133;
t55 = t238 + t254 + qJDD(5) + t158 * pkin(4) - qJ(5) * t185 + (t273 + t166) * qJD(3) + t275;
t268 = t193 * t55;
t267 = t196 * t55;
t40 = -t130 * pkin(5) - t143 * pkin(9) + t148 * t132 + t55;
t266 = t201 * t40;
t82 = t91 + t149;
t265 = t201 * t82;
t232 = t148 * t241;
t104 = t130 + t232;
t77 = t193 * t104 + t196 * t222;
t264 = t202 * t77;
t263 = t204 * t40;
t262 = t204 * t82;
t260 = t177 * t201;
t259 = t177 * t204;
t108 = t119 + t157;
t257 = t193 * t108;
t256 = t196 * t108;
t246 = -t144 - t184;
t161 = (t188 + t189) * qJDD(2);
t162 = t184 + t185;
t245 = pkin(2) * t162 + pkin(8) * t161;
t240 = qJD(6) + t177;
t236 = t202 * t91;
t235 = t202 * t119;
t9 = t201 * t11 + t204 * t12;
t88 = t202 * t96 - t127;
t89 = t205 * t96 + t254;
t44 = t202 * t88 + t205 * t89;
t227 = -t204 * t130 + t201 * t131;
t15 = t193 * t33 + t196 * t225;
t13 = t202 * t15 + t205 * t55;
t16 = -t193 * t225 + t196 * t33;
t219 = t201 * t130 + t204 * t131;
t218 = t205 * (-t184 + t207) + t253;
t215 = pkin(2) - t221;
t214 = t205 * t270 - pkin(2) - t255;
t213 = (-qJD(6) + t177) * t115 - t227;
t80 = -t113 * qJD(6) + t219;
t212 = qJD(3) * t273 + t275;
t73 = t228 * t202 + t226;
t210 = t212 + t238;
t74 = t229 * t241 + t211;
t163 = t184 - t185;
t141 = t202 * t157;
t135 = -t144 + t184;
t134 = t143 - t184;
t122 = t202 * t181 + t141;
t121 = (t158 - t231) * t205;
t118 = t205 * t156 + t202 * t159;
t117 = -t184 - t143;
t116 = (t161 * t203 + t162 * t206) * t195;
t105 = t131 - t233;
t103 = -t130 + t232;
t102 = t177 * t113;
t101 = -t143 - t144;
t98 = -t112 + t175;
t97 = t111 - t175;
t94 = -t112 - t175;
t93 = -t193 * t246 - t256;
t92 = t196 * t246 - t257;
t90 = t112 - t111;
t86 = -t175 - t111;
t85 = t196 * t117 - t285;
t84 = t193 * t117 + t284;
t79 = -t115 * qJD(6) - t227;
t78 = t196 * t104 - t193 * t222;
t76 = (-t113 * t204 + t115 * t201) * t177;
t75 = (-t113 * t201 - t115 * t204) * t177;
t72 = t210 + t254;
t71 = -t111 - t112;
t70 = t205 * t105 + t202 * t92;
t68 = t205 * t103 + t202 * t84;
t67 = t102 + t80;
t66 = -t102 + t80;
t65 = -t240 * t113 + t219;
t62 = t240 * t115 + t227;
t59 = t204 * t97 - t265;
t58 = -t201 * t98 + t282;
t57 = t201 * t97 + t262;
t56 = t204 * t98 + t283;
t52 = -t115 * t260 + t204 * t80;
t51 = t115 * t259 + t201 * t80;
t50 = t113 * t259 - t201 * t79;
t49 = t113 * t260 + t204 * t79;
t48 = t205 * t101 + t264;
t47 = -t201 * t94 - t262;
t46 = t204 * t94 - t265;
t43 = t204 * t86 - t283;
t42 = t201 * t86 + t282;
t39 = t202 * t73 + t205 * t72;
t38 = -t201 * t66 - t204 * t62;
t37 = t201 * t67 + t204 * t213;
t36 = -t201 * t62 + t204 * t66;
t35 = t201 * t213 - t204 * t67;
t30 = -t193 * t46 + t196 * t47;
t29 = t193 * t47 + t196 * t46;
t28 = -t193 * t42 + t196 * t43;
t27 = t193 * t43 + t196 * t42;
t25 = -pkin(9) * t46 + t263;
t23 = -pkin(9) * t42 + t266;
t22 = t202 * t29 + t205 * t65;
t21 = t202 * t27 + t205 * t62;
t20 = -pkin(5) * t65 + pkin(9) * t47 + t266;
t19 = -pkin(5) * t62 + pkin(9) * t43 - t263;
t18 = -t193 * t35 + t196 * t37;
t17 = t193 * t37 + t196 * t35;
t14 = t202 * t17 + t205 * t71;
t6 = -pkin(5) * t40 + pkin(9) * t9;
t5 = -pkin(9) * t35 - t8;
t4 = -pkin(5) * t71 + pkin(9) * t37 + t9;
t3 = t196 * t9 - t272;
t2 = t193 * t9 + t271;
t1 = t202 * t2 + t205 * t40;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t190, 0, 0, 0, 0, 0, 0, (qJDD(2) * t206 - t203 * t208) * t195, (-qJDD(2) * t203 - t206 * t208) * t195, 0, t198 * t133 + (t110 * t203 - t206 * t224) * t195, 0, 0, 0, 0, 0, 0, -t288, -t291, t116, t198 * (t202 * t89 - t205 * t88) + (t203 * t44 - t206 * t95) * t195, 0, 0, 0, 0, 0, 0, t116, t288, t291, t198 * (t202 * t72 - t205 * t73) + (t203 * t39 - t206 * t74) * t195, 0, 0, 0, 0, 0, 0, t198 * (t202 * t103 - t205 * t84) + (t203 * t68 - t206 * t85) * t195, t198 * (t202 * t105 - t205 * t92) + (t203 * t70 - t206 * t93) * t195, t198 * (t202 * t101 - t205 * t77) + (t203 * t48 - t206 * t78) * t195, t198 * (-t205 * t15 + t202 * t55) + (t203 * t13 - t206 * t16) * t195, 0, 0, 0, 0, 0, 0, t198 * (t202 * t62 - t205 * t27) + (t203 * t21 - t206 * t28) * t195, t198 * (t202 * t65 - t205 * t29) + (t203 * t22 - t206 * t30) * t195, t198 * (-t205 * t17 + t202 * t71) + (t203 * t14 - t206 * t18) * t195, t198 * (-t205 * t2 + t202 * t40) + (t203 * t1 - t206 * t3) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t224, -t110, 0, 0, t122, t118, t218, t121, t274, 0, pkin(2) * t159 - t205 * t95 - t286, -pkin(2) * t156 + t202 * t95 - t290, t44 + t245, -pkin(2) * t95 + pkin(8) * t44, 0, -t218, -t274, t122, t118, t121, (pkin(3) * t162 + t210) * t205 + (qJ(4) * t162 + t127 + t73) * t202 + t245, -t159 * t215 + t205 * t74 + t286, t202 * (-pkin(3) * t231 + t241 * t273 - t211) + t290 + t215 * t156, pkin(8) * t39 - t215 * t74, t235 + t205 * (-t193 * t131 - t196 * t232), t202 * (t144 - t143) + t205 * (t193 * t103 - t196 * t105), -t202 * t222 + t205 * (-t196 * t135 - t285), -t235 + t205 * (-t196 * t130 - t193 * t233), t202 * t104 + t205 * (-t193 * t134 - t256), t141 + t205 * (t146 * t193 + t148 * t196) * t241, t202 * (pkin(4) * t84 + t225) + t205 * (pkin(4) * t103 + t267) + pkin(8) * t68 + t214 * t85, t202 * (pkin(4) * t92 - t33) + t205 * (pkin(4) * t105 - t268) + pkin(8) * t70 + t214 * t93, pkin(4) * t264 + t205 * (pkin(4) * t101 - t16) + pkin(8) * t48 + t214 * t78, t16 * t214 + (pkin(4) + pkin(8)) * t13, t236 + t205 * (-t193 * t52 - t196 * t51), t202 * t90 + t205 * (-t193 * t38 - t196 * t36), t202 * t67 + t205 * (-t193 * t58 - t196 * t56), -t236 + t205 * (-t193 * t50 - t196 * t49), t202 * t213 + t205 * (-t193 * t59 - t196 * t57), t202 * t149 + t205 * (-t193 * t76 - t196 * t75), t202 * (pkin(4) * t27 + pkin(5) * t42 - t11) + t205 * (pkin(4) * t62 - t196 * t19 - t193 * t23) + pkin(8) * t21 + t214 * t28, t202 * (pkin(4) * t29 + pkin(5) * t46 - t12) + t205 * (pkin(4) * t65 - t193 * t25 - t196 * t20) + pkin(8) * t22 + t214 * t30, t202 * (pkin(4) * t17 + pkin(5) * t35) + t205 * (pkin(4) * t71 - t193 * t5 - t196 * t4) + pkin(8) * t14 + t214 * t18, t202 * (pkin(4) * t2 + pkin(5) * t8) + t205 * (pkin(4) * t40 + pkin(9) * t272 - t196 * t6) + pkin(8) * t1 + t214 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t163, t183, t234, t237, qJDD(3), -t88, -t89, 0, 0, qJDD(3), -t183, -t237, -t234, t163, t234, (-pkin(3) * t202 + qJ(4) * t205) * qJDD(2), -pkin(3) * t167 - qJ(4) * t174 + t73, -pkin(3) * t172 + t254 + (qJDD(3) + t168) * qJ(4) + t212, -pkin(3) * t73 + qJ(4) * t72, t196 * t131 - t193 * t232, -t196 * t103 - t193 * t105, -t193 * t135 + t284, -t193 * t130 + t196 * t233, t196 * t134 - t257, (-t146 * t196 + t148 * t193) * t241, qJ(4) * t103 + t270 * t84 + t268, qJ(4) * t105 + t270 * t92 + t267, qJ(4) * t101 + t270 * t77 - t15, qJ(4) * t55 + t270 * t15, -t193 * t51 + t196 * t52, -t193 * t36 + t196 * t38, -t193 * t56 + t196 * t58, -t193 * t49 + t196 * t50, -t193 * t57 + t196 * t59, -t193 * t75 + t196 * t76, qJ(4) * t62 - t193 * t19 + t196 * t23 + t270 * t27, qJ(4) * t65 - t193 * t20 + t196 * t25 + t270 * t29, qJ(4) * t71 + t270 * t17 - t193 * t4 + t196 * t5, -pkin(9) * t271 + qJ(4) * t40 - t193 * t6 + t270 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t167, t172, t73, 0, 0, 0, 0, 0, 0, t84, t92, t77, t15, 0, 0, 0, 0, 0, 0, t27, t29, t17, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t105, t101, t55, 0, 0, 0, 0, 0, 0, t62, t65, t71, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90, t67, -t91, t213, t149, -t11, -t12, 0, 0;];
tauJ_reg  = t7;
