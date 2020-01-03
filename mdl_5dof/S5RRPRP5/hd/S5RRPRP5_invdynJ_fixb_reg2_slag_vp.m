% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:55:02
% EndTime: 2019-12-31 19:55:07
% DurationCPUTime: 3.02s
% Computational Cost: add. (4678->370), mult. (11221->446), div. (0->0), fcn. (8267->12), ass. (0->191)
t165 = sin(qJ(4));
t162 = sin(pkin(8));
t163 = cos(pkin(8));
t168 = cos(qJ(2));
t231 = t163 * t168;
t217 = qJD(1) * t231;
t166 = sin(qJ(2));
t226 = qJD(1) * t166;
t106 = -t162 * t226 + t217;
t259 = cos(qJ(4));
t114 = t162 * t168 + t163 * t166;
t265 = t114 * qJD(1);
t195 = -t165 * t106 - t259 * t265;
t222 = t168 * qJDD(1);
t223 = t166 * qJDD(1);
t199 = t162 * t223 - t163 * t222;
t270 = t265 * qJD(2);
t70 = t199 + t270;
t224 = qJD(1) * qJD(2);
t215 = t166 * t224;
t186 = t114 * qJDD(1) - t162 * t215;
t214 = t168 * t224;
t71 = t163 * t214 + t186;
t176 = t195 * qJD(4) - t165 * t71 - t259 * t70;
t158 = qJD(2) + qJD(4);
t271 = t158 * t195;
t278 = t176 - t271;
t216 = qJD(4) * t259;
t225 = qJD(4) * t165;
t194 = t106 * t216 - t165 * t70 - t225 * t265 + t259 * t71;
t60 = t259 * t106 - t165 * t265;
t239 = t60 * t158;
t277 = t194 - t239;
t250 = t60 ^ 2;
t276 = t195 ^ 2;
t202 = t276 - t250;
t252 = t168 * pkin(2);
t146 = pkin(1) + t252;
t123 = -t146 * qJD(1) + qJD(3);
t77 = -t106 * pkin(3) + t123;
t27 = -pkin(4) * t60 + qJ(5) * t195 + t77;
t275 = t27 * t60;
t274 = t77 * t60;
t273 = t195 * t77;
t272 = t60 * t195;
t156 = qJDD(2) + qJDD(4);
t151 = t156 * pkin(4);
t264 = qJDD(5) - t151;
t269 = -t195 * t27 + t264;
t34 = -pkin(4) * t195 - t60 * qJ(5);
t145 = t163 * pkin(2) + pkin(3);
t255 = t106 * pkin(7);
t246 = qJ(3) + pkin(6);
t132 = t246 * t168;
t120 = qJD(1) * t132;
t111 = t163 * t120;
t131 = t246 * t166;
t119 = qJD(1) * t131;
t75 = t162 * t119 - t111;
t189 = t75 - t255;
t258 = pkin(2) * t162;
t221 = t165 * t258;
t254 = t265 * pkin(7);
t109 = t162 * t120;
t76 = -t163 * t119 - t109;
t49 = t76 - t254;
t244 = -qJD(4) * t221 + t145 * t216 - t165 * t189 - t259 * t49;
t159 = qJ(2) + pkin(8);
t154 = qJ(4) + t159;
t143 = sin(t154);
t144 = cos(t154);
t268 = t144 * pkin(4) + t143 * qJ(5);
t101 = t165 * t145 + t259 * t258;
t147 = t156 * qJ(5);
t150 = t158 * qJD(5);
t267 = t147 + t150;
t167 = sin(qJ(1));
t169 = cos(qJ(1));
t266 = g(1) * t167 - g(2) * t169;
t78 = -t163 * t131 - t162 * t132;
t53 = -t114 * pkin(7) + t78;
t197 = t162 * t166 - t231;
t79 = -t162 * t131 + t163 * t132;
t54 = -t197 * pkin(7) + t79;
t33 = t165 * t53 + t259 * t54;
t190 = t197 * qJD(2);
t209 = qJD(2) * t246;
t102 = t168 * qJD(3) - t166 * t209;
t103 = -t166 * qJD(3) - t168 * t209;
t51 = -t162 * t102 + t163 * t103;
t174 = pkin(7) * t190 + t51;
t192 = t114 * qJD(2);
t52 = t163 * t102 + t162 * t103;
t42 = -pkin(7) * t192 + t52;
t8 = t165 * t174 + t53 * t216 - t54 * t225 + t259 * t42;
t263 = t143 * t266 + t33 * t156 + t8 * t158;
t261 = t265 ^ 2;
t170 = qJD(2) ^ 2;
t260 = t70 * pkin(3);
t256 = g(3) * t168;
t253 = t166 * pkin(2);
t67 = qJDD(2) * pkin(2) + t103 * qJD(1) - qJDD(1) * t131;
t74 = t102 * qJD(1) + qJDD(1) * t132;
t36 = t162 * t67 + t163 * t74;
t245 = qJD(5) + t244;
t243 = t101 * qJD(4) - t165 * t49 + t259 * t189;
t242 = qJD(2) * pkin(2);
t113 = -t119 + t242;
t68 = t163 * t113 - t109;
t45 = qJD(2) * pkin(3) - t254 + t68;
t69 = t162 * t113 + t111;
t48 = t69 + t255;
t25 = t165 * t45 + t259 * t48;
t240 = t25 * t158;
t237 = pkin(6) * qJDD(1);
t236 = t265 * t106;
t235 = t143 * t167;
t234 = t143 * t169;
t233 = t144 * t167;
t232 = t144 * t169;
t24 = -t165 * t48 + t259 * t45;
t230 = qJD(5) - t24;
t160 = t166 ^ 2;
t161 = t168 ^ 2;
t228 = t160 - t161;
t227 = t160 + t161;
t149 = t166 * t242;
t171 = qJD(1) ^ 2;
t219 = t166 * t171 * t168;
t218 = t243 * t195;
t153 = cos(t159);
t213 = pkin(3) * t153 + t252;
t80 = pkin(2) * t226 + pkin(3) * t265;
t152 = sin(t159);
t121 = -pkin(3) * t152 - t253;
t212 = -pkin(4) * t143 + t121;
t35 = -t162 * t74 + t163 * t67;
t23 = qJDD(2) * pkin(3) - t71 * pkin(7) + t35;
t26 = -t70 * pkin(7) + t36;
t3 = t165 * t23 + t45 * t216 - t48 * t225 + t259 * t26;
t208 = t165 * t26 + t48 * t216 + t45 * t225 - t259 * t23;
t207 = t166 * t214;
t206 = g(1) * t169 + g(2) * t167;
t73 = t259 * t114 - t165 * t197;
t38 = t73 * qJD(4) - t165 * t190 + t259 * t192;
t187 = t259 * t197;
t72 = t165 * t114 + t187;
t204 = -t176 * t72 - t38 * t60;
t203 = -t276 - t250;
t198 = t156 * t72 + t158 * t38;
t196 = -0.2e1 * pkin(1) * t224 - pkin(6) * qJDD(2);
t193 = -g(1) * t232 - g(2) * t233 - g(3) * t143 + t3;
t100 = t259 * t145 - t221;
t97 = pkin(2) * t215 - t146 * qJDD(1) + qJDD(3);
t188 = g(1) * t234 + g(2) * t235 - g(3) * t144 - t208;
t185 = t194 + t239;
t184 = t24 * t158 - t193;
t32 = t165 * t54 - t259 * t53;
t9 = t33 * qJD(4) + t165 * t42 - t259 * t174;
t183 = g(1) * t233 - g(2) * t232 - t32 * t156 - t9 * t158;
t81 = pkin(3) * t192 + t149;
t37 = t114 * t225 + t158 * t187 + t165 * t192;
t182 = t176 * t73 - t194 * t72 + t195 * t38 - t37 * t60;
t181 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t170 + t266;
t180 = pkin(1) * t171 + t206 - t237;
t46 = t97 + t260;
t84 = t197 * pkin(3) - t146;
t179 = t176 * t33 + t194 * t32 - t195 * t9 + t60 * t8 - t206;
t178 = -t266 + t97;
t177 = -t243 * t158 + t188;
t175 = -t188 + t269;
t173 = -t176 - t271;
t5 = -pkin(4) * t176 - qJ(5) * t194 + qJD(5) * t195 + t46;
t157 = -pkin(7) - t246;
t125 = qJ(5) * t232;
t124 = qJ(5) * t233;
t118 = pkin(1) + t213;
t112 = t169 * t118;
t104 = t106 ^ 2;
t98 = -pkin(4) - t100;
t96 = qJ(5) + t101;
t31 = t72 * pkin(4) - t73 * qJ(5) + t84;
t30 = t34 + t80;
t18 = t73 * t156 - t37 * t158;
t15 = t158 * qJ(5) + t25;
t14 = -t158 * pkin(4) + t230;
t10 = t38 * pkin(4) + t37 * qJ(5) - t73 * qJD(5) + t81;
t6 = t194 * t73 + t195 * t37;
t2 = t208 + t264;
t1 = t3 + t267;
t4 = [0, 0, 0, 0, 0, qJDD(1), t266, t206, 0, 0, t160 * qJDD(1) + 0.2e1 * t207, 0.2e1 * t166 * t222 - 0.2e1 * t228 * t224, qJDD(2) * t166 + t170 * t168, t161 * qJDD(1) - 0.2e1 * t207, qJDD(2) * t168 - t170 * t166, 0, t166 * t196 + t168 * t181, -t166 * t181 + t168 * t196, 0.2e1 * t227 * t237 - t206, -g(1) * (-t167 * pkin(1) + t169 * pkin(6)) - g(2) * (t169 * pkin(1) + t167 * pkin(6)) + (t227 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t71 * t114 - t190 * t265, -t114 * t70 - t71 * t197 + (-t197 * t106 - t114 * t265) * qJD(2), t114 * qJDD(2) - t197 * t170, -t106 * t192 + t70 * t197, -t197 * qJDD(2) - t114 * t170, 0, -t146 * t70 + t97 * t197 + t78 * qJDD(2) + t266 * t153 + (-t106 * t253 + t123 * t114 + t51) * qJD(2), -t79 * qJDD(2) + t97 * t114 - t146 * t71 - t266 * t152 + (-t123 * t197 + t253 * t265 - t52) * qJD(2), t52 * t106 - t79 * t70 - t36 * t197 - t51 * t265 - t78 * t71 - t35 * t114 + (-t114 * t69 + t197 * t68) * qJD(2) - t206, t36 * t79 + t69 * t52 + t35 * t78 + t68 * t51 - t97 * t146 + t123 * t149 - g(1) * (-t167 * t146 + t169 * t246) - g(2) * (t169 * t146 + t167 * t246), t6, t182, t18, t204, -t198, 0, -t176 * t84 + t77 * t38 + t46 * t72 - t60 * t81 + t183, t194 * t84 - t195 * t81 - t77 * t37 + t46 * t73 - t263, t208 * t73 + t24 * t37 - t25 * t38 - t3 * t72 + t179, t3 * t33 + t25 * t8 + t208 * t32 - t24 * t9 + t46 * t84 + t77 * t81 - g(1) * (-t167 * t118 - t169 * t157) - g(2) * (-t167 * t157 + t112), t6, t18, -t182, 0, t198, t204, -t10 * t60 - t176 * t31 + t27 * t38 + t5 * t72 + t183, -t1 * t72 - t14 * t37 - t15 * t38 + t2 * t73 + t179, t10 * t195 - t194 * t31 + t27 * t37 - t5 * t73 + t263, -g(2) * t112 + t1 * t33 + t27 * t10 + t14 * t9 + t15 * t8 + t2 * t32 + t5 * t31 + (g(1) * t157 - g(2) * t268) * t169 + (-g(1) * (-t118 - t268) + g(2) * t157) * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t228 * t171, t223, t219, t222, qJDD(2), t166 * t180 - t256, g(3) * t166 + t168 * t180, 0, 0, -t236, -t104 + t261, (-t106 + t217) * qJD(2) + t186, t236, -t199, qJDD(2), -g(3) * t153 - t75 * qJD(2) - t123 * t265 + t206 * t152 + (t163 * qJDD(2) + t106 * t226) * pkin(2) + t35, g(3) * t152 + t76 * qJD(2) - t123 * t106 + t206 * t153 + (-t162 * qJDD(2) - t226 * t265) * pkin(2) - t36, (t69 + t75) * t265 + (t68 - t76) * t106 + (-t70 * t162 - t71 * t163) * pkin(2), -t68 * t75 - t69 * t76 + (-t256 + t36 * t162 + t163 * t35 + (-qJD(1) * t123 + t206) * t166) * pkin(2), t272, t202, t277, -t272, t278, t156, t100 * t156 + t60 * t80 + t177 + t273, -t101 * t156 - t244 * t158 + t195 * t80 - t193 - t274, -t100 * t194 + t101 * t176 - t195 * t25 - t218 + (t24 + t244) * t60, -g(3) * t213 - t100 * t208 + t3 * t101 - t206 * t121 - t243 * t24 + t244 * t25 - t77 * t80, t272, t277, -t202, t156, -t278, -t272, -t98 * t156 + t30 * t60 + t177 - t269, -t15 * t195 + t194 * t98 + t96 * t176 - t218 + (-t14 + t245) * t60, t96 * t156 + t245 * t158 - t195 * t30 + t193 + t267 + t275, t1 * t96 + t2 * t98 - t27 * t30 - g(1) * (t169 * t212 + t125) - g(2) * (t167 * t212 + t124) - g(3) * (t213 + t268) + t245 * t15 + t243 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199 + 0.2e1 * t270, (t106 + t217) * qJD(2) + t186, -t104 - t261, -t69 * t106 + t265 * t68 + t178, 0, 0, 0, 0, 0, 0, t173, t185, t203, -t195 * t24 - t25 * t60 + t178 + t260, 0, 0, 0, 0, 0, 0, t173, t203, -t185, t14 * t195 - t15 * t60 - t266 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t202, t277, -t272, t278, t156, t188 + t240 + t273, t184 - t274, 0, 0, t272, t277, -t202, t156, -t278, -t272, t34 * t60 + t151 - t175 + t240, -pkin(4) * t194 + t176 * qJ(5) - (t15 - t25) * t195 - (t14 - t230) * t60, -t195 * t34 + 0.2e1 * t147 + 0.2e1 * t150 - t184 + t275, t1 * qJ(5) - t2 * pkin(4) - t27 * t34 - t14 * t25 - g(1) * (-pkin(4) * t234 + t125) - g(2) * (-pkin(4) * t235 + t124) - g(3) * t268 + t230 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156 + t272, t277, -t158 ^ 2 - t276, -t15 * t158 + t175;];
tau_reg = t4;
