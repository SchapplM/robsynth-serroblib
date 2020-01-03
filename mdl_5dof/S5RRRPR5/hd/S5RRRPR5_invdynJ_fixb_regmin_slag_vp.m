% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:20
% EndTime: 2019-12-31 21:14:26
% DurationCPUTime: 2.49s
% Computational Cost: add. (4382->326), mult. (10656->441), div. (0->0), fcn. (7951->14), ass. (0->197)
t166 = qJ(2) + qJ(3);
t158 = sin(t166);
t159 = cos(t166);
t172 = sin(qJ(1));
t176 = cos(qJ(1));
t207 = g(1) * t176 + g(2) * t172;
t284 = -g(3) * t159 + t158 * t207;
t157 = pkin(9) + t166;
t147 = cos(t157);
t161 = qJDD(2) + qJDD(3);
t167 = sin(pkin(9));
t168 = cos(pkin(9));
t170 = sin(qJ(3));
t174 = cos(qJ(3));
t175 = cos(qJ(2));
t238 = qJD(1) * t175;
t171 = sin(qJ(2));
t239 = qJD(1) * t171;
t101 = -t170 * t238 - t174 * t239;
t272 = pkin(6) + pkin(7);
t127 = t272 * t175;
t120 = qJD(1) * t127;
t106 = t174 * t120;
t126 = t272 * t171;
t118 = qJD(1) * t126;
t255 = qJD(2) * pkin(2);
t108 = -t118 + t255;
t200 = -t108 * t170 - t106;
t233 = qJD(1) * qJD(2);
t222 = t175 * t233;
t232 = t171 * qJDD(1);
t280 = t222 + t232;
t79 = qJDD(2) * pkin(2) - t272 * t280;
t223 = t171 * t233;
t231 = t175 * qJDD(1);
t81 = t272 * (-t223 + t231);
t183 = qJD(3) * t200 - t170 * t81 + t174 * t79;
t163 = qJD(2) + qJD(3);
t225 = t174 * t238;
t226 = t170 * t239;
t59 = qJD(3) * t225 - t163 * t226 + t170 * t231 + t174 * t280;
t18 = pkin(3) * t161 - qJ(4) * t59 + qJD(4) * t101 + t183;
t237 = qJD(3) * t170;
t274 = (qJD(3) * t108 + t81) * t174 - t120 * t237 + t170 * t79;
t203 = t170 * t232 - t174 * t231;
t113 = t170 * t175 + t171 * t174;
t78 = t163 * t113;
t60 = qJD(1) * t78 + t203;
t99 = -t225 + t226;
t22 = -qJ(4) * t60 - qJD(4) * t99 + t274;
t5 = -t167 * t22 + t168 * t18;
t3 = -pkin(4) * t161 - t5;
t283 = g(3) * t147 + t3;
t216 = t101 * t167 - t168 * t99;
t275 = qJD(5) - t216;
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t201 = -t101 * t168 - t167 * t99;
t61 = -t163 * t173 + t169 * t201;
t282 = t275 * t61;
t213 = t173 * t275;
t34 = -t167 * t59 - t168 * t60;
t31 = qJDD(5) - t34;
t279 = -t169 * t31 - t275 * t213;
t160 = t175 * pkin(2);
t260 = pkin(1) + t160;
t148 = pkin(3) * t167 + pkin(8);
t271 = pkin(3) * t101;
t42 = pkin(4) * t201 - pkin(8) * t216 - t271;
t278 = (qJD(5) * t148 + t42) * t275;
t155 = pkin(2) * t239;
t153 = pkin(2) * t174 + pkin(3);
t249 = t168 * t170;
t94 = pkin(2) * t249 + t153 * t167;
t89 = pkin(8) + t94;
t277 = (qJD(5) * t89 + t155 + t42) * t275;
t102 = t170 * t120;
t215 = t108 * t174 - t102;
t95 = t101 * qJ(4);
t57 = t215 + t95;
t242 = -t126 * t170 + t127 * t174;
t227 = qJD(2) * t272;
t119 = t171 * t227;
t121 = t175 * t227;
t182 = -qJD(3) * t242 + t119 * t170 - t121 * t174;
t112 = t170 * t171 - t174 * t175;
t77 = t163 * t112;
t180 = qJ(4) * t77 - qJD(4) * t113 + t182;
t236 = qJD(3) * t174;
t189 = -t119 * t174 - t121 * t170 - t126 * t236 - t127 * t237;
t38 = -qJ(4) * t78 - qJD(4) * t112 + t189;
t15 = t167 * t180 + t168 * t38;
t125 = t260 * qJD(1);
t80 = pkin(3) * t99 + qJD(4) - t125;
t37 = -pkin(4) * t216 - pkin(8) * t201 + t80;
t6 = t167 * t18 + t168 * t22;
t221 = pkin(8) * t161 + qJD(5) * t37 + t6;
t256 = qJ(4) * t99;
t58 = -t200 - t256;
t253 = t167 * t58;
t52 = pkin(3) * t163 + t57;
t27 = t168 * t52 - t253;
t25 = -pkin(4) * t163 - t27;
t208 = pkin(3) * t112 - t260;
t74 = t112 * t168 + t113 * t167;
t75 = -t112 * t167 + t113 * t168;
t43 = pkin(4) * t74 - pkin(8) * t75 + t208;
t209 = -t126 * t174 - t127 * t170;
t193 = -qJ(4) * t113 + t209;
t68 = -qJ(4) * t112 + t242;
t45 = t167 * t193 + t168 * t68;
t50 = -t167 * t78 - t168 * t77;
t273 = -(qJD(5) * t43 + t15) * t275 - t221 * t74 + t25 * t50 + t3 * t75 - t45 * t31;
t266 = t25 * t216;
t265 = t25 * t75;
t264 = t43 * t31;
t263 = t61 * t201;
t63 = t163 * t169 + t173 * t201;
t262 = t63 * t201;
t261 = t275 * t201;
t53 = t168 * t58;
t28 = t167 * t52 + t53;
t210 = t118 * t170 - t106;
t196 = t210 + t256;
t257 = pkin(2) * qJD(3);
t243 = -t118 * t174 - t102;
t64 = t95 + t243;
t259 = -t167 * t64 + t168 * t196 + (t167 * t174 + t249) * t257;
t250 = t167 * t170;
t258 = -t167 * t196 - t168 * t64 + (t168 * t174 - t250) * t257;
t254 = t101 * t99;
t234 = qJD(5) * t173;
t235 = qJD(5) * t169;
t35 = -t167 * t60 + t168 * t59;
t20 = t161 * t169 + t163 * t234 + t173 * t35 - t201 * t235;
t252 = t169 * t20;
t248 = t169 * t172;
t247 = t169 * t176;
t246 = t172 * t173;
t245 = t173 * t176;
t241 = pkin(3) * t159 + t160;
t164 = t171 ^ 2;
t240 = -t175 ^ 2 + t164;
t156 = t171 * t255;
t230 = t75 * t235;
t224 = pkin(3) * t78 + t156;
t96 = pkin(2) * t223 - qJDD(1) * t260;
t184 = t60 * pkin(3) + qJDD(4) + t96;
t11 = -pkin(4) * t34 - pkin(8) * t35 + t184;
t26 = pkin(8) * t163 + t28;
t219 = qJD(5) * t26 - t11;
t214 = -t161 * t173 + t169 * t35;
t206 = g(1) * t172 - g(2) * t176;
t205 = t201 * t28 + t216 * t27;
t204 = t275 * t50 + t31 * t75;
t13 = t169 * t37 + t173 * t26;
t202 = t13 * t201 + t169 * t283 + t25 * t234;
t12 = -t169 * t26 + t173 * t37;
t146 = sin(t157);
t199 = -t12 * t201 + t25 * t235 + (g(1) * t245 + g(2) * t246) * t146;
t198 = t173 * t31 + (t169 * t216 - t235) * t275;
t197 = g(3) * t146 - t221;
t93 = -pkin(2) * t250 + t153 * t168;
t195 = t207 * t146;
t194 = -0.2e1 * pkin(1) * t233 - pkin(6) * qJDD(2);
t33 = t168 * t57 - t253;
t188 = -t148 * t31 + t275 * t33 - t266;
t187 = -t258 * t275 - t89 * t31 - t266;
t177 = qJD(2) ^ 2;
t186 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t177 + t206;
t178 = qJD(1) ^ 2;
t185 = pkin(1) * t178 - pkin(6) * qJDD(1) + t207;
t181 = g(3) * t158 - t125 * t99 + t159 * t207 - t274;
t179 = -t125 * t101 + t183 + t284;
t162 = -qJ(4) - t272;
t149 = -pkin(3) * t168 - pkin(4);
t117 = pkin(1) + t241;
t88 = -pkin(4) - t93;
t87 = t147 * t245 + t248;
t86 = -t147 * t247 + t246;
t85 = -t147 * t246 + t247;
t84 = t147 * t248 + t245;
t65 = t101 ^ 2 - t99 ^ 2;
t49 = -t167 * t77 + t168 * t78;
t48 = -t203 + (-qJD(1) * t113 - t101) * t163;
t47 = t163 * t99 + t59;
t44 = t167 * t68 - t168 * t193;
t32 = t167 * t57 + t53;
t21 = qJD(5) * t63 + t214;
t17 = pkin(4) * t49 - pkin(8) * t50 + t224;
t14 = t167 * t38 - t168 * t180;
t10 = t173 * t11;
t9 = t213 * t63 + t252;
t8 = -t262 - t279;
t7 = t198 + t263;
t1 = (t20 - t282) * t173 + (-t275 * t63 - t21) * t169;
t2 = [qJDD(1), t206, t207, qJDD(1) * t164 + 0.2e1 * t171 * t222, 0.2e1 * t171 * t231 - 0.2e1 * t233 * t240, qJDD(2) * t171 + t175 * t177, qJDD(2) * t175 - t171 * t177, 0, t171 * t194 + t175 * t186, -t171 * t186 + t175 * t194, t101 * t77 + t113 * t59, t101 * t78 - t112 * t59 - t113 * t60 + t77 * t99, t113 * t161 - t163 * t77, -t112 * t161 - t163 * t78, 0, t96 * t112 - t125 * t78 + t156 * t99 + t159 * t206 + t161 * t209 + t163 * t182 - t260 * t60, -t101 * t156 + t96 * t113 + t125 * t77 - t158 * t206 - t161 * t242 - t163 * t189 - t260 * t59, t14 * t201 + t15 * t216 - t27 * t50 - t28 * t49 + t34 * t45 + t35 * t44 - t5 * t75 - t6 * t74 - t207, t6 * t45 + t28 * t15 - t5 * t44 - t27 * t14 + t184 * t208 + t80 * t224 - g(1) * (-t117 * t172 - t162 * t176) - g(2) * (t117 * t176 - t162 * t172), -t63 * t230 + (t20 * t75 + t50 * t63) * t173, (-t169 * t63 - t173 * t61) * t50 + (-t252 - t173 * t21 + (t169 * t61 - t173 * t63) * qJD(5)) * t75, t173 * t204 + t20 * t74 - t230 * t275 + t49 * t63, -t234 * t275 * t75 - t169 * t204 - t21 * t74 - t49 * t61, t275 * t49 + t31 * t74, -g(1) * t85 - g(2) * t87 + t10 * t74 + t12 * t49 + t14 * t61 + t44 * t21 + (t17 * t275 + t264 + (-t26 * t74 - t275 * t45 + t265) * qJD(5)) * t173 + t273 * t169, -g(1) * t84 - g(2) * t86 - t13 * t49 + t14 * t63 + t44 * t20 + (-(-qJD(5) * t45 + t17) * t275 - t264 + t219 * t74 - qJD(5) * t265) * t169 + t273 * t173; 0, 0, 0, -t171 * t178 * t175, t240 * t178, t232, t231, qJDD(2), -g(3) * t175 + t171 * t185, g(3) * t171 + t175 * t185, -t254, t65, t47, t48, t161, -t210 * t163 + (t161 * t174 - t163 * t237 - t239 * t99) * pkin(2) + t179, t243 * t163 + (t101 * t239 - t161 * t170 - t163 * t236) * pkin(2) + t181, t201 * t259 + t216 * t258 + t34 * t94 - t35 * t93 + t205, t6 * t94 + t5 * t93 - t80 * (t155 - t271) - g(3) * t241 + t258 * t28 - t259 * t27 - t207 * (-pkin(2) * t171 - pkin(3) * t158), t9, t1, t8, t7, -t261, t88 * t21 + t259 * t61 + (-t283 - t277) * t173 + t187 * t169 + t199, t88 * t20 + t259 * t63 + t187 * t173 + (-t195 + t277) * t169 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, t65, t47, t48, t161, -t163 * t200 + t179, t163 * t215 + t181, -t32 * t201 - t33 * t216 + (t167 * t34 - t168 * t35) * pkin(3) + t205, t27 * t32 - t28 * t33 + (t101 * t80 + t167 * t6 + t168 * t5 + t284) * pkin(3), t9, t1, t8, t7, -t261, t149 * t21 - t32 * t61 + t188 * t169 + (-t283 - t278) * t173 + t199, t149 * t20 - t32 * t63 + t188 * t173 + (-t195 + t278) * t169 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201 ^ 2 - t216 ^ 2, t201 * t27 - t216 * t28 + t184 - t206, 0, 0, 0, 0, 0, t198 - t263, -t262 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t61, -t61 ^ 2 + t63 ^ 2, t20 + t282, -t214 + (-qJD(5) + t275) * t63, t31, -g(1) * t86 + g(2) * t84 + t13 * t275 + t169 * t197 - t234 * t26 - t25 * t63 + t10, g(1) * t87 - g(2) * t85 + t12 * t275 + t169 * t219 + t173 * t197 + t25 * t61;];
tau_reg = t2;
