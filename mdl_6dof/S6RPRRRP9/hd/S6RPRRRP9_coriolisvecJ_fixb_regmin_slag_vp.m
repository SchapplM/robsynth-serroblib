% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:31
% EndTime: 2019-03-09 06:28:40
% DurationCPUTime: 3.13s
% Computational Cost: add. (4083->356), mult. (8881->503), div. (0->0), fcn. (5809->6), ass. (0->183)
t154 = sin(qJ(4));
t155 = sin(qJ(3));
t205 = t155 * qJD(1);
t194 = t154 * t205;
t250 = pkin(8) + pkin(9);
t195 = qJD(4) * t250;
t158 = cos(qJ(3));
t175 = pkin(3) * t158 + pkin(8) * t155;
t119 = t175 * qJD(1);
t157 = cos(qJ(4));
t159 = -pkin(1) - pkin(7);
t138 = t159 * qJD(1) + qJD(2);
t224 = t158 * t138;
t237 = t154 * t119 + t157 * t224;
t261 = pkin(9) * t194 + t154 * t195 + t237;
t101 = t157 * t119;
t197 = t154 * t224;
t227 = t155 * t157;
t199 = pkin(9) * t227;
t265 = t157 * t195 - t197 + t101 + (pkin(4) * t158 + t199) * qJD(1);
t204 = t157 * qJD(3);
t215 = qJD(1) * t158;
t111 = -t154 * t215 + t204;
t188 = t157 * t215;
t206 = t154 * qJD(3);
t112 = t188 + t206;
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t172 = t153 * t111 + t156 * t112;
t251 = t172 ^ 2;
t61 = -t156 * t111 + t153 * t112;
t59 = t61 ^ 2;
t264 = -t59 + t251;
t201 = qJD(4) + qJD(5);
t207 = qJD(5) * t156;
t210 = qJD(4) * t157;
t226 = t156 * t157;
t230 = t153 * t154;
t241 = t153 * t194 - t156 * t210 - t157 * t207 + t201 * t230 - t205 * t226;
t115 = t153 * t157 + t156 * t154;
t68 = t201 * t115;
t94 = t115 * qJD(1);
t240 = t155 * t94 + t68;
t263 = t172 * t61;
t262 = t61 * qJ(6);
t209 = qJD(4) * t158;
t190 = t157 * t209;
t193 = t155 * t206;
t260 = t190 - t193;
t143 = qJD(4) + t205;
t136 = qJD(5) + t143;
t202 = qJD(3) * qJD(4);
t219 = qJD(4) * t188 + t154 * t202;
t167 = qJD(1) * t193 - t219;
t208 = qJD(5) * t153;
t189 = t155 * t204;
t191 = t154 * t209;
t166 = -t189 - t191;
t74 = t166 * qJD(1) + t157 * t202;
t23 = -t111 * t207 + t112 * t208 - t153 * t167 - t156 * t74;
t259 = t61 * t136 - t23;
t122 = t155 * t138;
t102 = qJD(3) * pkin(8) + t122;
t174 = t155 * pkin(3) - t158 * pkin(8);
t123 = qJ(2) + t174;
t96 = t123 * qJD(1);
t239 = t154 * t96;
t55 = t157 * t102 + t239;
t109 = t175 * qJD(3) + qJD(2);
t83 = t109 * qJD(1);
t77 = t157 * t83;
t164 = -t55 * qJD(4) + t77;
t213 = qJD(3) * t158;
t17 = -t74 * pkin(9) + (pkin(4) * qJD(1) - t138 * t154) * t213 + t164;
t192 = t158 * t204;
t212 = qJD(4) * t154;
t169 = -t102 * t212 + t138 * t192 + t154 * t83 + t96 * t210;
t21 = t167 * pkin(9) + t169;
t54 = -t154 * t102 + t157 * t96;
t43 = -t112 * pkin(9) + t54;
t38 = t143 * pkin(4) + t43;
t44 = t111 * pkin(9) + t55;
t181 = -t153 * t17 - t156 * t21 - t38 * t207 + t44 * t208;
t103 = -qJD(3) * pkin(3) - t224;
t70 = -t111 * pkin(4) + t103;
t258 = t70 * t61 + t181;
t200 = 0.2e1 * qJD(1);
t31 = t61 * pkin(5) + qJD(6) + t70;
t257 = t172 * t31;
t256 = t265 * t156;
t255 = qJ(6) * t172;
t114 = -t226 + t230;
t86 = t114 * t155;
t84 = t115 * t155;
t131 = t250 * t154;
t132 = t250 * t157;
t220 = -t153 * t131 + t156 * t132;
t152 = t158 ^ 2;
t217 = t155 ^ 2 - t152;
t254 = t131 * t207 + t132 * t208 + t265 * t153 + t261 * t156;
t42 = t156 * t44;
t12 = t153 * t38 + t42;
t186 = -t153 * t21 + t156 * t17;
t163 = -t12 * qJD(5) + t186;
t253 = -t70 * t172 + t163;
t24 = t172 * qJD(5) + t153 * t74 - t156 * t167;
t252 = t136 * t172 - t24;
t40 = t153 * t44;
t11 = t156 * t38 - t40;
t6 = t11 - t255;
t5 = t136 * pkin(5) + t6;
t249 = t5 - t6;
t248 = -qJ(6) * t240 - t114 * qJD(6) - t254;
t247 = -pkin(5) * t215 + qJ(6) * t241 - t220 * qJD(5) - t115 * qJD(6) + t153 * t261 - t256;
t87 = t114 * t158;
t246 = -qJD(3) * t87 - t201 * t84 - t94;
t245 = t114 * qJD(1) - t115 * t213 + t201 * t86;
t244 = t156 * t43 - t40;
t108 = t157 * t123;
t228 = t154 * t159;
t187 = pkin(4) - t228;
t225 = t157 * t158;
t58 = -pkin(9) * t225 + t187 * t155 + t108;
t135 = t159 * t227;
t218 = t154 * t123 + t135;
t229 = t154 * t158;
t69 = -pkin(9) * t229 + t218;
t242 = t153 * t58 + t156 * t69;
t238 = t74 * t154;
t236 = t103 * t154;
t235 = t103 * t155;
t234 = t111 * t143;
t233 = t143 * t154;
t232 = t143 * t155;
t231 = t143 * t157;
t160 = qJD(3) ^ 2;
t223 = t160 * t155;
t222 = t160 * t158;
t161 = qJD(1) ^ 2;
t221 = t161 * qJ(2);
t216 = -t160 - t161;
t214 = qJD(3) * t155;
t211 = qJD(4) * t155;
t203 = qJD(1) * qJD(3);
t198 = qJD(2) * t200;
t196 = t143 * t229;
t150 = -t157 * pkin(4) - pkin(3);
t145 = t158 * t203;
t90 = t157 * t109;
t28 = t90 + (-t135 + (pkin(9) * t158 - t123) * t154) * qJD(4) + (t187 * t158 + t199) * qJD(3);
t165 = t154 * t109 + t123 * t210 + t159 * t192 - t211 * t228;
t30 = -pkin(9) * t260 + t165;
t185 = -t153 * t30 + t156 * t28;
t184 = -t153 * t43 - t42;
t183 = -t153 * t69 + t156 * t58;
t180 = -t156 * t131 - t153 * t132;
t110 = pkin(4) * t229 - t158 * t159;
t179 = -t111 + t204;
t178 = -t112 + t206;
t177 = qJD(1) + t211;
t176 = -t122 + (t194 + t212) * pkin(4);
t171 = t178 * t155;
t170 = t153 * t28 + t156 * t30 + t58 * t207 - t69 * t208;
t75 = pkin(4) * t260 + t159 * t214;
t51 = -t167 * pkin(4) + t138 * t214;
t10 = t24 * pkin(5) + t51;
t149 = t156 * pkin(4) + pkin(5);
t133 = t155 * t145;
t85 = t115 * t158;
t50 = -t114 * qJ(6) + t220;
t49 = -t115 * qJ(6) + t180;
t37 = -t208 * t229 + (t201 * t225 - t193) * t156 + t166 * t153;
t35 = -t153 * t193 + t156 * t189 + t68 * t158;
t25 = -t85 * qJ(6) + t242;
t22 = t155 * pkin(5) + t87 * qJ(6) + t183;
t9 = t244 - t255;
t8 = t184 + t262;
t7 = t12 - t262;
t4 = -t37 * qJ(6) - t85 * qJD(6) + t170;
t3 = pkin(5) * t213 + t35 * qJ(6) - qJD(5) * t242 + t87 * qJD(6) + t185;
t2 = -t24 * qJ(6) - t61 * qJD(6) - t181;
t1 = pkin(5) * t145 + t23 * qJ(6) - qJD(6) * t172 + t163;
t13 = [0, 0, 0, 0, t198, qJ(2) * t198, -0.2e1 * t133, 0.2e1 * t217 * t203, -t223, -t222, 0, -t159 * t223 + (qJ(2) * t213 + qJD(2) * t155) * t200, -t159 * t222 + (-qJ(2) * t214 + qJD(2) * t158) * t200, t166 * t112 + t74 * t225 (-t157 * t111 + t112 * t154) * t214 + (t157 * t167 - t238 + (-t154 * t111 - t112 * t157) * qJD(4)) * t158, -t143 * t191 + t74 * t155 + (t112 * t158 + (qJD(1) * t152 - t232) * t157) * qJD(3), -t143 * t190 - t219 * t155 + (t111 * t158 + (t217 * qJD(1) + t232) * t154) * qJD(3), t143 * t213 + t133 (-t123 * t212 + t90) * t143 + (-t159 * t219 + t103 * t210 + (t108 * qJD(1) - t143 * t228 + t54) * qJD(3)) * t158 + (t77 + (-t111 * t159 - t236) * qJD(3) + (-t239 + (-t143 * t159 - t102) * t157) * qJD(4)) * t155, -t165 * t143 - t169 * t155 + (-t103 * t212 - t159 * t74) * t158 + ((-qJD(1) * t218 - t55) * t158 + (t159 * t112 + (-t103 + t224) * t157) * t155) * qJD(3), -t172 * t35 + t23 * t87, -t172 * t37 + t23 * t85 + t87 * t24 + t35 * t61, -t35 * t136 - t23 * t155 + (-qJD(1) * t87 + t172) * t213, -t37 * t136 - t24 * t155 + (-qJD(1) * t85 - t61) * t213, t136 * t213 + t133, t185 * t136 + t186 * t155 + t75 * t61 + t110 * t24 + t51 * t85 + t70 * t37 + (-t12 * t155 - t136 * t242) * qJD(5) + (qJD(1) * t183 + t11) * t213, -t170 * t136 + t181 * t155 + t75 * t172 - t110 * t23 - t51 * t87 - t70 * t35 + (-t242 * qJD(1) - t12) * t213, t1 * t87 - t172 * t3 - t2 * t85 + t22 * t23 - t25 * t24 + t5 * t35 - t7 * t37 - t4 * t61, t2 * t25 + t7 * t4 + t1 * t22 + t5 * t3 + t10 * (t85 * pkin(5) + t110) + t31 * (t37 * pkin(5) + t75); 0, 0, 0, 0, -t161, -t221, 0, 0, 0, 0, 0, t216 * t155, t216 * t158, 0, 0, 0, 0, 0, -t158 * t219 - t177 * t231 + (-t155 * t111 - t196) * qJD(3), -t158 * t74 + t177 * t233 + (-t143 * t225 + (t112 - t188) * t155) * qJD(3), 0, 0, 0, 0, 0, -t158 * t24 + t245 * t136 + (t155 * t61 - t215 * t84) * qJD(3), t158 * t23 - t246 * t136 + (t155 * t172 + t215 * t86) * qJD(3), -t172 * t245 - t84 * t23 + t86 * t24 - t246 * t61, -t1 * t84 - t10 * t158 - t2 * t86 + t31 * t214 + t245 * t5 + t246 * t7; 0, 0, 0, 0, 0, 0, t158 * t161 * t155, -t217 * t161, 0, 0, 0, -t158 * t221, t155 * t221, t112 * t231 + t238 (t74 + t234) * t157 + (qJD(1) * t171 - t112 * qJD(4) - t219) * t154, t143 * t210 + (t143 * t227 + t158 * t178) * qJD(1), -t143 * t212 + (-t154 * t232 + t158 * t179) * qJD(1), -t143 * t215, -pkin(3) * t219 - t101 * t143 + (-t155 * t179 + t196) * t138 + (-pkin(8) * t231 + t236) * qJD(4) + (-t54 * t158 + (qJD(3) * t174 + t235) * t154) * qJD(1), -pkin(3) * t74 + t237 * t143 + t138 * t171 + (pkin(8) * t233 + t103 * t157) * qJD(4) + (t55 * t158 + (-pkin(8) * t213 + t235) * t157) * qJD(1), -t23 * t115 - t172 * t241, t23 * t114 - t115 * t24 - t172 * t240 + t241 * t61, -t241 * t136 + (qJD(3) * t115 - t172) * t215, -t240 * t136 + (-qJD(3) * t114 + t61) * t215, -t136 * t215, t51 * t114 + t150 * t24 + t240 * t70 + t176 * t61 + (-t132 * t207 + (qJD(5) * t131 + t261) * t153 - t256) * t136 + (qJD(3) * t180 - t11) * t215, t51 * t115 - t150 * t23 - t241 * t70 + t176 * t172 + t254 * t136 + (-qJD(3) * t220 + t12) * t215, -t1 * t115 - t2 * t114 - t172 * t247 + t49 * t23 - t50 * t24 - t240 * t7 + t241 * t5 - t248 * t61, t2 * t50 + t1 * t49 + t10 * (t114 * pkin(5) + t150) + t248 * t7 + t247 * t5 + (pkin(4) * t233 + t240 * pkin(5) - t122) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 * t111, -t111 ^ 2 + t112 ^ 2, t74 - t234, t112 * t143 + t167, t145, -qJD(3) * t197 - t103 * t112 + t55 * t143 + t164, -t103 * t111 + t54 * t143 - t169, t263, t264, t259, t252, t145, -t184 * t136 + (-t112 * t61 - t136 * t208 + t145 * t156) * pkin(4) + t253, t244 * t136 + (-t112 * t172 - t136 * t207 - t145 * t153) * pkin(4) + t258, t149 * t23 - t5 * t61 + t7 * t172 + t9 * t61 + t8 * t172 + (-t153 * t24 + (t153 * t172 - t156 * t61) * qJD(5)) * pkin(4), -pkin(5) * t257 + t1 * t149 - t5 * t8 - t7 * t9 + (-t31 * t112 + t2 * t153 + (-t153 * t5 + t156 * t7) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, t264, t259, t252, t145, t12 * t136 + t253, t11 * t136 + t258, pkin(5) * t23 - t249 * t61, t249 * t7 + (t1 - t257) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59 - t251, t172 * t5 + t7 * t61 + t10;];
tauc_reg  = t13;
