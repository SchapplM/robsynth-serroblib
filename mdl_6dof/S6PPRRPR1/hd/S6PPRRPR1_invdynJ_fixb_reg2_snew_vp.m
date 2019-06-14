% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:05:25
% EndTime: 2019-05-04 20:05:34
% DurationCPUTime: 4.02s
% Computational Cost: add. (20863->365), mult. (40209->573), div. (0->0), fcn. (32369->16), ass. (0->235)
t188 = sin(pkin(7));
t192 = cos(pkin(7));
t187 = sin(pkin(12));
t191 = cos(pkin(12));
t189 = sin(pkin(6));
t193 = cos(pkin(6));
t255 = sin(pkin(11));
t256 = cos(pkin(11));
t209 = t255 * g(1) - t256 * g(2);
t241 = -g(3) + qJDD(1);
t207 = t189 * t241 + t193 * t209;
t210 = -t256 * g(1) - t255 * g(2);
t204 = -t187 * t210 + t191 * t207;
t206 = -t189 * t209 + t193 * t241 + qJDD(2);
t273 = t188 * t206 + t192 * t204;
t186 = sin(pkin(13));
t190 = cos(pkin(13));
t196 = sin(qJ(4));
t240 = qJD(3) * t196;
t156 = -t190 * qJD(4) + t186 * t240;
t158 = t186 * qJD(4) + t190 * t240;
t195 = sin(qJ(6));
t198 = cos(qJ(6));
t129 = t198 * t156 + t195 * t158;
t199 = cos(qJ(4));
t239 = t199 * qJD(3);
t173 = -qJD(6) + t239;
t118 = t129 * t173;
t237 = qJD(3) * qJD(4);
t232 = t199 * t237;
t236 = t196 * qJDD(3);
t162 = t232 + t236;
t141 = t186 * qJDD(4) + t190 * t162;
t230 = -t190 * qJDD(4) + t186 * t162;
t99 = -t129 * qJD(6) + t198 * t141 - t195 * t230;
t272 = t118 + t99;
t176 = t196 * t237;
t235 = t199 * qJDD(3);
t163 = -t176 + t235;
t253 = t158 * t156;
t211 = -t163 - t253;
t271 = t186 * t211;
t270 = t190 * t211;
t159 = -qJDD(6) + t163;
t131 = -t195 * t156 + t198 * t158;
t254 = t131 * t129;
t208 = -t159 - t254;
t269 = t195 * t208;
t268 = t198 * t208;
t146 = t156 * t239;
t123 = -t141 + t146;
t147 = t158 * t239;
t121 = -t230 - t147;
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t119 = t187 * t207 + t191 * t210;
t228 = t197 * t119 - t200 * t273;
t90 = t200 * t119 + t197 * t273;
t267 = t197 * t90 - t200 * t228;
t231 = t195 * t141 + t198 * t230;
t80 = (qJD(6) + t173) * t131 + t231;
t127 = t129 ^ 2;
t128 = t131 ^ 2;
t266 = t156 ^ 2;
t155 = t158 ^ 2;
t171 = t173 ^ 2;
t265 = qJD(4) ^ 2;
t227 = -t199 * pkin(4) - t196 * qJ(5);
t160 = t227 * qJD(3);
t202 = -t188 * t204 + t192 * t206;
t201 = qJD(3) ^ 2;
t88 = -t201 * pkin(3) + qJDD(3) * pkin(9) + t90;
t61 = t196 * t202 + t199 * t88;
t58 = -t265 * pkin(4) + qJDD(4) * qJ(5) + t160 * t239 + t61;
t226 = t162 + t232;
t87 = -qJDD(3) * pkin(3) - t201 * pkin(9) + t228;
t64 = -t226 * qJ(5) + (-t163 + t176) * pkin(4) + t87;
t39 = 0.2e1 * qJD(5) * t158 + t186 * t58 - t190 * t64;
t30 = t211 * pkin(5) + t123 * pkin(10) - t39;
t142 = -pkin(5) * t239 - t158 * pkin(10);
t40 = -0.2e1 * qJD(5) * t156 + t186 * t64 + t190 * t58;
t31 = -t266 * pkin(5) - t230 * pkin(10) + t142 * t239 + t40;
t18 = t195 * t31 - t198 * t30;
t19 = t195 * t30 + t198 * t31;
t10 = -t198 * t18 + t195 * t19;
t264 = t186 * t10;
t108 = t199 * t202;
t57 = qJDD(5) - t108 - t265 * qJ(5) - qJDD(4) * pkin(4) + (qJD(3) * t160 + t88) * t196;
t263 = t186 * t57;
t262 = t190 * t10;
t261 = t190 * t57;
t49 = t230 * pkin(5) - t266 * pkin(10) + t158 * t142 + t57;
t260 = t195 * t49;
t258 = t198 * t49;
t252 = t173 * t195;
t251 = t173 * t198;
t124 = t163 - t253;
t250 = t186 * t124;
t249 = t189 * t187;
t248 = t189 * t191;
t247 = t190 * t124;
t246 = t191 * t192;
t100 = t159 - t254;
t245 = t195 * t100;
t172 = t196 * t201 * t199;
t167 = qJDD(4) + t172;
t244 = t196 * t167;
t243 = t198 * t100;
t168 = qJDD(4) - t172;
t242 = t199 * t168;
t234 = t199 * t254;
t233 = t199 * t253;
t11 = t195 * t18 + t198 * t19;
t22 = t186 * t39 + t190 * t40;
t60 = t196 * t88 - t108;
t42 = t196 * t60 + t199 * t61;
t5 = t190 * t11 - t264;
t3 = t196 * t49 + t199 * t5;
t4 = t186 * t11 + t262;
t229 = t197 * t3 - t200 * t4;
t16 = t196 * t57 + t199 * t22;
t21 = t186 * t40 - t190 * t39;
t225 = t16 * t197 - t200 * t21;
t84 = -t118 + t99;
t53 = -t195 * t80 - t198 * t84;
t55 = t195 * t84 - t198 * t80;
t28 = -t186 * t53 + t190 * t55;
t93 = -t127 - t128;
t24 = t196 * t93 + t199 * t28;
t27 = t186 * t55 + t190 * t53;
t224 = t197 * t24 - t200 * t27;
t107 = -t171 - t127;
t65 = t195 * t107 + t268;
t66 = t198 * t107 - t269;
t46 = -t186 * t65 + t190 * t66;
t79 = (qJD(6) - t173) * t131 + t231;
t33 = t196 * t79 + t199 * t46;
t45 = t186 * t66 + t190 * t65;
t223 = t197 * t33 - t200 * t45;
t112 = -t128 - t171;
t67 = t198 * t112 + t245;
t68 = -t195 * t112 + t243;
t51 = -t186 * t67 + t190 * t68;
t37 = t196 * t272 + t199 * t51;
t50 = t186 * t68 + t190 * t67;
t222 = t197 * t37 - t200 * t50;
t221 = t197 * t42 - t200 * t87;
t116 = -t155 - t266;
t97 = t190 * t121 - t186 * t123;
t70 = t196 * t116 + t199 * t97;
t96 = t186 * t121 + t190 * t123;
t220 = t197 * t70 - t200 * t96;
t183 = t199 ^ 2;
t180 = t183 * t201;
t133 = -t180 - t266;
t103 = t186 * t133 + t270;
t104 = t190 * t133 - t271;
t120 = -t147 + t230;
t86 = t199 * t104 + t196 * t120;
t219 = -t103 * t200 + t197 * t86;
t145 = -t155 - t180;
t110 = t190 * t145 + t250;
t111 = -t186 * t145 + t247;
t122 = t141 + t146;
t92 = t199 * t111 + t196 * t122;
t218 = -t110 * t200 + t197 * t92;
t170 = -t180 - t265;
t138 = t199 * t170 - t244;
t164 = -0.2e1 * t176 + t235;
t217 = t138 * t197 + t164 * t200;
t182 = t196 ^ 2;
t179 = t182 * t201;
t169 = -t179 - t265;
t139 = -t196 * t169 - t242;
t161 = 0.2e1 * t232 + t236;
t216 = t139 * t197 - t161 * t200;
t165 = (t182 + t183) * qJDD(3);
t166 = t179 + t180;
t215 = t165 * t197 + t166 * t200;
t214 = t200 * qJDD(3) - t197 * t201;
t213 = -t197 * qJDD(3) - t200 * t201;
t212 = -pkin(3) + t227;
t153 = t199 * t163;
t150 = t214 * t188;
t149 = t213 * t188;
t144 = -t155 + t180;
t143 = -t180 + t266;
t137 = -t196 * t168 + t199 * t169;
t136 = t199 * t167 + t196 * t170;
t132 = t215 * t188;
t115 = -t128 + t171;
t114 = t127 - t171;
t109 = t128 - t127;
t106 = t192 * t137 + t188 * t216;
t105 = t192 * t136 + t188 * t217;
t98 = -t131 * qJD(6) - t231;
t95 = (t129 * t198 - t131 * t195) * t173;
t94 = (t129 * t195 + t131 * t198) * t173;
t91 = t196 * t111 - t199 * t122;
t85 = t196 * t104 - t199 * t120;
t78 = t198 * t114 + t245;
t77 = -t195 * t115 + t268;
t76 = t195 * t114 - t243;
t75 = t198 * t115 + t269;
t74 = t131 * t252 + t198 * t99;
t73 = -t131 * t251 + t195 * t99;
t72 = -t129 * t251 - t195 * t98;
t71 = -t129 * t252 + t198 * t98;
t69 = -t199 * t116 + t196 * t97;
t54 = -t195 * t272 - t198 * t79;
t52 = -t195 * t79 + t198 * t272;
t48 = t188 * t218 + t192 * t91;
t47 = t267 * t188 + t192 * t202;
t44 = t188 * t219 + t192 * t85;
t43 = t188 * t220 + t192 * t69;
t41 = t196 * t61 - t199 * t60;
t36 = t196 * t51 - t199 * t272;
t35 = -pkin(10) * t67 + t258;
t34 = -pkin(10) * t65 + t260;
t32 = t196 * t46 - t199 * t79;
t26 = -pkin(5) * t272 + pkin(10) * t68 + t260;
t25 = -pkin(5) * t79 + pkin(10) * t66 - t258;
t23 = t196 * t28 - t199 * t93;
t20 = t188 * t221 + t192 * t41;
t15 = t196 * t22 - t199 * t57;
t14 = t188 * t222 + t192 * t36;
t13 = t188 * t223 + t192 * t32;
t12 = t188 * t224 + t192 * t23;
t9 = -pkin(10) * t53 - t10;
t8 = -pkin(5) * t49 + pkin(10) * t11;
t7 = -pkin(5) * t93 + pkin(10) * t55 + t11;
t6 = t192 * t15 + t188 * t225;
t2 = t196 * t5 - t199 * t49;
t1 = t229 * t188 + t192 * t2;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t249 + t193 * t206 + t204 * t248, 0, 0, 0, 0, 0, 0, t193 * t150 + (t187 * t213 + t214 * t246) * t189, t193 * t149 + (-t187 * t214 + t213 * t246) * t189, 0, (t197 * t228 + t200 * t90) * t249 + (-t188 * t202 + t267 * t192) * t248 + t193 * t47, 0, 0, 0, 0, 0, 0, t193 * t105 + (t187 * (t200 * t138 - t197 * t164) + t191 * (-t188 * t136 + t192 * t217)) * t189, t193 * t106 + (t187 * (t200 * t139 + t197 * t161) + t191 * (-t188 * t137 + t192 * t216)) * t189, t193 * t132 + (t187 * (t200 * t165 - t197 * t166) + t215 * t246) * t189, t193 * t20 + (t187 * (t197 * t87 + t200 * t42) + t191 * (-t188 * t41 + t192 * t221)) * t189, 0, 0, 0, 0, 0, 0, t193 * t44 + (t187 * (t197 * t103 + t200 * t86) + t191 * (-t188 * t85 + t192 * t219)) * t189, t193 * t48 + (t187 * (t197 * t110 + t200 * t92) + t191 * (-t188 * t91 + t192 * t218)) * t189, t193 * t43 + (t187 * (t197 * t96 + t200 * t70) + t191 * (-t188 * t69 + t192 * t220)) * t189, t193 * t6 + (t187 * (t200 * t16 + t197 * t21) + t191 * (-t188 * t15 + t192 * t225)) * t189, 0, 0, 0, 0, 0, 0, t193 * t13 + (t187 * (t197 * t45 + t200 * t33) + t191 * (-t188 * t32 + t192 * t223)) * t189, t193 * t14 + (t187 * (t197 * t50 + t200 * t37) + t191 * (-t188 * t36 + t192 * t222)) * t189, t193 * t12 + (t187 * (t197 * t27 + t200 * t24) + t191 * (-t188 * t23 + t192 * t224)) * t189, t193 * t1 + (t187 * (t197 * t4 + t200 * t3) + t191 * (-t188 * t2 + t192 * t229)) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, 0, 0, 0, 0, 0, t150, t149, 0, t47, 0, 0, 0, 0, 0, 0, t105, t106, t132, t20, 0, 0, 0, 0, 0, 0, t44, t48, t43, t6, 0, 0, 0, 0, 0, 0, t13, t14, t12, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t228, -t90, 0, 0, t226 * t196, t199 * t161 + t196 * t164, t244 + t199 * (-t179 + t265), -t196 * t232 + t153, t196 * (t180 - t265) + t242, 0, pkin(3) * t164 + pkin(9) * t138 - t199 * t87, -pkin(3) * t161 + pkin(9) * t139 + t196 * t87, pkin(3) * t166 + pkin(9) * t165 + t42, -pkin(3) * t87 + pkin(9) * t42, t196 * (t190 * t141 + t186 * t147) - t233, t196 * (-t190 * t120 - t186 * t122) + t199 * (-t155 + t266), t196 * (-t186 * t144 + t270) + t199 * t123, t196 * (-t190 * t146 + t186 * t230) + t233, t196 * (t190 * t143 + t250) - t199 * t121, t153 + t196 * (t156 * t190 - t158 * t186) * t239, t196 * (-qJ(5) * t103 + t263) + t199 * (-pkin(4) * t103 + t39) - pkin(3) * t103 + pkin(9) * t86, t196 * (-qJ(5) * t110 + t261) + t199 * (-pkin(4) * t110 + t40) - pkin(3) * t110 + pkin(9) * t92, pkin(9) * t70 - t196 * t21 + t212 * t96, pkin(9) * t16 + t21 * t212, t196 * (-t186 * t73 + t190 * t74) - t234, t196 * (-t186 * t52 + t190 * t54) - t199 * t109, t196 * (-t186 * t75 + t190 * t77) - t199 * t84, t196 * (-t186 * t71 + t190 * t72) + t234, t196 * (-t186 * t76 + t190 * t78) + t199 * t80, t196 * (-t186 * t94 + t190 * t95) + t199 * t159, t196 * (-qJ(5) * t45 - t186 * t25 + t190 * t34) + t199 * (-pkin(4) * t45 - pkin(5) * t65 + t18) - pkin(3) * t45 + pkin(9) * t33, t196 * (-qJ(5) * t50 - t186 * t26 + t190 * t35) + t199 * (-pkin(4) * t50 - pkin(5) * t67 + t19) - pkin(3) * t50 + pkin(9) * t37, t196 * (-qJ(5) * t27 - t186 * t7 + t190 * t9) + t199 * (-pkin(4) * t27 - pkin(5) * t53) - pkin(3) * t27 + pkin(9) * t24, t196 * (-pkin(10) * t262 - qJ(5) * t4 - t186 * t8) + t199 * (-pkin(4) * t4 - pkin(5) * t10) - pkin(3) * t4 + pkin(9) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t179 - t180, t236, t172, t235, qJDD(4), -t60, -t61, 0, 0, t186 * t141 - t190 * t147, -t186 * t120 + t190 * t122, t190 * t144 + t271, -t186 * t146 - t190 * t230, t186 * t143 - t247, (t156 * t186 + t158 * t190) * t239, -pkin(4) * t120 + qJ(5) * t104 - t261, -pkin(4) * t122 + qJ(5) * t111 + t263, -pkin(4) * t116 + qJ(5) * t97 + t22, -pkin(4) * t57 + qJ(5) * t22, t186 * t74 + t190 * t73, t186 * t54 + t190 * t52, t186 * t77 + t190 * t75, t186 * t72 + t190 * t71, t186 * t78 + t190 * t76, t186 * t95 + t190 * t94, -pkin(4) * t79 + qJ(5) * t46 + t186 * t34 + t190 * t25, -pkin(4) * t272 + qJ(5) * t51 + t186 * t35 + t190 * t26, -pkin(4) * t93 + qJ(5) * t28 + t186 * t9 + t190 * t7, -pkin(4) * t49 - pkin(10) * t264 + qJ(5) * t5 + t190 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t122, t116, t57, 0, 0, 0, 0, 0, 0, t79, t272, t93, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, t109, t84, -t254, -t80, -t159, -t18, -t19, 0, 0;];
tauJ_reg  = t17;
