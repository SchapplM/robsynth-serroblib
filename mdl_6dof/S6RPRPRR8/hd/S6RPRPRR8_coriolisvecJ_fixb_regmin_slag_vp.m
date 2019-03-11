% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:31
% EndTime: 2019-03-09 03:59:40
% DurationCPUTime: 2.72s
% Computational Cost: add. (3697->335), mult. (8341->471), div. (0->0), fcn. (6130->8), ass. (0->181)
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t219 = sin(pkin(10));
t220 = cos(pkin(10));
t162 = t220 * t149 + t219 * t152;
t242 = t162 * qJD(1);
t251 = qJD(5) + t242;
t115 = -t219 * t149 + t220 * t152;
t110 = t115 * qJD(1);
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t195 = t151 * qJD(3);
t84 = t148 * t110 - t195;
t255 = t251 * t84;
t147 = sin(qJ(6));
t150 = cos(qJ(6));
t86 = t148 * qJD(3) + t151 * t110;
t169 = t147 * t84 - t150 * t86;
t30 = t147 * t86 + t150 * t84;
t254 = t169 * t30;
t201 = qJD(5) * t148;
t215 = t148 * t242;
t253 = t201 + t215;
t177 = t251 * t151;
t98 = t110 * qJD(3);
t229 = t148 * t98;
t252 = -t177 * t251 - t229;
t250 = t169 ^ 2 - t30 ^ 2;
t198 = qJD(6) * t150;
t199 = qJD(6) * t147;
t99 = qJD(3) * t242;
t43 = qJD(5) * t195 - t110 * t201 - t151 * t99;
t228 = t148 * t99;
t44 = t86 * qJD(5) - t228;
t8 = -t147 * t44 + t150 * t43 - t84 * t198 - t86 * t199;
t97 = qJD(6) + t251;
t249 = t30 * t97 + t8;
t153 = -pkin(1) - pkin(7);
t125 = t153 * qJD(1) + qJD(2);
t205 = qJD(1) * t149;
t102 = -qJ(4) * t205 + t149 * t125;
t182 = t220 * t102;
t204 = qJD(1) * t152;
t103 = -qJ(4) * t204 + t152 * t125;
t94 = qJD(3) * pkin(3) + t103;
t51 = t219 * t94 + t182;
t42 = qJD(3) * pkin(8) + t51;
t120 = pkin(3) * t205 + qJD(1) * qJ(2) + qJD(4);
t52 = pkin(4) * t242 - t110 * pkin(8) + t120;
t21 = t148 * t52 + t151 * t42;
t14 = -t84 * pkin(9) + t21;
t11 = t14 * t199;
t91 = t219 * t102;
t50 = t220 * t94 - t91;
t41 = -qJD(3) * pkin(4) - t50;
t24 = t84 * pkin(5) + t41;
t248 = t24 * t30 + t11;
t157 = t169 * qJD(6) - t147 * t43 - t150 * t44;
t247 = -t169 * t97 + t157;
t194 = t152 * qJD(4);
t203 = qJD(3) * t149;
t156 = -t125 * t203 + (qJ(4) * t203 - t194) * qJD(1);
t196 = t149 * qJD(4);
t202 = qJD(3) * t152;
t83 = t125 * t202 + (-qJ(4) * t202 - t196) * qJD(1);
t28 = t219 * t156 + t220 * t83;
t143 = qJD(1) * qJD(2);
t191 = qJD(1) * qJD(3);
t186 = t152 * t191;
t208 = pkin(3) * t186 + t143;
t39 = t98 * pkin(4) + t99 * pkin(8) + t208;
t35 = t151 * t39;
t158 = -t21 * qJD(5) - t148 * t28 + t35;
t2 = t98 * pkin(5) - t43 * pkin(9) + t158;
t200 = qJD(5) * t151;
t167 = t148 * t39 + t151 * t28 + t52 * t200 - t42 * t201;
t5 = -t44 * pkin(9) + t167;
t189 = -t147 * t5 + t150 * t2;
t20 = -t148 * t42 + t151 * t52;
t13 = -t86 * pkin(9) + t20;
t10 = pkin(5) * t251 + t13;
t227 = t150 * t14;
t4 = t147 * t10 + t227;
t246 = -t4 * qJD(6) + t24 * t169 + t189;
t117 = t147 * t148 - t150 * t151;
t118 = t147 * t151 + t150 * t148;
t241 = qJD(5) + qJD(6);
t77 = t241 * t118;
t245 = -t117 * t98 - t77 * t97;
t66 = t118 * t115;
t178 = qJD(3) * t219;
t179 = qJD(3) * t220;
t109 = -t149 * t179 - t152 * t178;
t214 = t151 * t109;
t243 = -t115 * t201 + t214;
t230 = t118 * t98;
t76 = t241 * t117;
t234 = -t117 * t242 - t76;
t240 = -t234 * t97 - t230;
t239 = 0.2e1 * t143;
t209 = qJ(4) - t153;
t121 = t209 * t149;
t122 = t209 * t152;
t79 = -t220 * t121 - t219 * t122;
t238 = t79 * t98;
t134 = t219 * pkin(3) + pkin(8);
t237 = pkin(9) + t134;
t62 = t220 * t103 - t91;
t64 = pkin(3) * t204 + t110 * pkin(4) + pkin(8) * t242;
t236 = t148 * t64 + t151 * t62;
t210 = t149 * pkin(3) + qJ(2);
t71 = pkin(4) * t162 - t115 * pkin(8) + t210;
t72 = t151 * t79;
t235 = t148 * t71 + t72;
t59 = t118 * t242;
t233 = t77 + t59;
t232 = t110 * t30;
t231 = t110 * t84;
t73 = t162 * t98;
t89 = t151 * t98;
t27 = -t220 * t156 + t219 * t83;
t226 = t27 * t115;
t225 = t27 * t151;
t224 = t169 * t110;
t223 = t43 * t148;
t222 = t86 * t110;
t108 = t149 * t178 - t152 * t179;
t221 = t97 * t108;
t218 = t109 * t148;
t217 = t115 * t148;
t216 = t115 * t151;
t154 = qJD(3) ^ 2;
t213 = t154 * t149;
t212 = t154 * t152;
t155 = qJD(1) ^ 2;
t211 = t155 * qJ(2);
t207 = t149 ^ 2 - t152 ^ 2;
t206 = -t154 - t155;
t197 = t120 * qJD(1);
t193 = pkin(3) * t202 + qJD(2);
t190 = 0.2e1 * qJD(1);
t188 = t115 * t200;
t185 = qJD(6) * t10 + t5;
t183 = qJD(5) * t237;
t176 = -qJD(5) * t162 - qJD(1);
t175 = -t59 * t97 + t245;
t61 = t219 * t103 + t182;
t174 = t253 * pkin(5) - t61;
t135 = -t220 * pkin(3) - pkin(4);
t113 = t237 * t148;
t171 = pkin(9) * t215 + qJD(6) * t113 + t148 * t183 + t236;
t114 = t237 * t151;
t58 = t151 * t64;
t170 = t110 * pkin(5) + qJD(6) * t114 - t148 * t62 + t58 + (pkin(9) * t242 + t183) * t151;
t100 = t209 * t203 - t194;
t101 = -qJD(3) * t122 - t196;
t55 = -t220 * t100 + t219 * t101;
t78 = -t219 * t121 + t220 * t122;
t168 = -t251 * t253 + t89;
t56 = t219 * t100 + t220 * t101;
t63 = -t108 * pkin(4) - t109 * pkin(8) + t193;
t166 = t148 * t63 + t151 * t56 + t71 * t200 - t79 * t201;
t165 = t188 + t218;
t164 = t97 * t118;
t163 = -t134 * t98 + t251 * t41;
t160 = -t51 * t108 + t50 * t109 + t162 * t28 - t226;
t123 = -t151 * pkin(5) + t135;
t69 = t151 * t71;
t67 = t117 * t115;
t54 = t151 * t63;
t48 = pkin(5) * t217 + t78;
t23 = t165 * pkin(5) + t55;
t22 = -pkin(9) * t217 + t235;
t19 = pkin(5) * t162 - pkin(9) * t216 - t148 * t79 + t69;
t17 = t44 * pkin(5) + t27;
t16 = -t199 * t217 + (t241 * t216 + t218) * t150 + t243 * t147;
t15 = -t117 * t109 - t241 * t66;
t7 = -t165 * pkin(9) + t166;
t6 = -pkin(9) * t214 - t108 * pkin(5) - t148 * t56 + t54 + (-t72 + (pkin(9) * t115 - t71) * t148) * qJD(5);
t3 = t150 * t10 - t147 * t14;
t1 = [0, 0, 0, 0, t239, qJ(2) * t239, -0.2e1 * t149 * t186, 0.2e1 * t207 * t191, -t213, -t212, 0, -t153 * t213 + (qJ(2) * t202 + qJD(2) * t149) * t190, -t153 * t212 + (-qJ(2) * t203 + qJD(2) * t152) * t190, t55 * t110 - t242 * t56 - t78 * t99 - t160 - t238, t120 * t193 + t208 * t210 + t27 * t78 + t28 * t79 - t50 * t55 + t51 * t56, t86 * t214 + (t43 * t151 - t86 * t201) * t115 (-t148 * t86 - t151 * t84) * t109 + (-t223 - t151 * t44 + (t148 * t84 - t151 * t86) * qJD(5)) * t115, -t86 * t108 + t162 * t43 + t98 * t216 + t243 * t251, t84 * t108 - t162 * t44 - t165 * t251 - t98 * t217, -t108 * t251 + t73 (-t200 * t79 + t54) * t251 + t69 * t98 + (-t200 * t42 + t35) * t162 - t20 * t108 + t55 * t84 + t78 * t44 + t41 * t188 + ((-qJD(5) * t71 - t56) * t251 - t238 + (-qJD(5) * t52 - t28) * t162 + t226 + t41 * t109) * t148, -t166 * t251 - t235 * t98 - t167 * t162 + t21 * t108 + t55 * t86 + t78 * t43 + t41 * t214 + (-t201 * t41 + t225) * t115, -t15 * t169 - t8 * t67, -t15 * t30 - t157 * t67 + t16 * t169 - t8 * t66, t108 * t169 + t15 * t97 + t162 * t8 - t67 * t98, t30 * t108 + t157 * t162 - t16 * t97 - t66 * t98, t73 - t221 (-t147 * t7 + t150 * t6) * t97 + (-t147 * t22 + t150 * t19) * t98 + t189 * t162 - t3 * t108 + t23 * t30 - t48 * t157 + t17 * t66 + t24 * t16 + ((-t147 * t19 - t150 * t22) * t97 - t4 * t162) * qJD(6), t4 * t108 + t11 * t162 + t24 * t15 - t17 * t67 - t23 * t169 + t48 * t8 + (-(-qJD(6) * t22 + t6) * t97 - t19 * t98 - t2 * t162) * t147 + (-(qJD(6) * t19 + t7) * t97 - t22 * t98 - t185 * t162) * t150; 0, 0, 0, 0, -t155, -t211, 0, 0, 0, 0, 0, t206 * t149, t206 * t152, t108 * t242 - t109 * t110 + t115 * t99 - t73, t160 - t197, 0, 0, 0, 0, 0, -t162 * t229 - t109 * t84 - t115 * t44 + (t108 * t148 + t151 * t176) * t251, -t162 * t89 - t109 * t86 - t115 * t43 + (t108 * t151 - t148 * t176) * t251, 0, 0, 0, 0, 0, -t109 * t30 + t115 * t157 + t108 * t164 + t117 * t97 * qJD(1) - (-t76 * t97 + t230) * t162, qJD(1) * t164 + t109 * t169 - t115 * t8 - t117 * t221 - t162 * t245; 0, 0, 0, 0, 0, 0, t152 * t155 * t149, -t207 * t155, 0, 0, 0, -t152 * t211, t149 * t211 (t51 - t61) * t110 - (t50 - t62) * t242 + (-t219 * t98 + t220 * t99) * pkin(3), t50 * t61 - t51 * t62 + (-t152 * t197 + t219 * t28 - t220 * t27) * pkin(3), t86 * t177 + t223 (t43 - t255) * t151 + (-t251 * t86 - t44) * t148, -t222 - t252, t168 + t231, -t251 * t110, -t20 * t110 + t135 * t44 - t225 - t61 * t84 + (-t134 * t200 - t58) * t251 + (t251 * t62 + t163) * t148, t21 * t110 + t135 * t43 + t27 * t148 - t61 * t86 + (t134 * t201 + t236) * t251 + t163 * t151, t8 * t118 - t169 * t234, -t8 * t117 + t118 * t157 + t169 * t233 - t234 * t30, t224 - t240, t175 + t232, -t97 * t110 (-t150 * t113 - t147 * t114) * t98 - t123 * t157 + t17 * t117 - t3 * t110 + (t147 * t171 - t150 * t170) * t97 + t174 * t30 + t233 * t24 -(-t147 * t113 + t150 * t114) * t98 + t123 * t8 + t17 * t118 + t4 * t110 + (t147 * t170 + t150 * t171) * t97 - t174 * t169 + t234 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 ^ 2 - t242 ^ 2, t50 * t110 + t242 * t51 + t208, 0, 0, 0, 0, 0, t168 - t231, -t222 + t252, 0, 0, 0, 0, 0, t175 - t232, t224 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t84, -t84 ^ 2 + t86 ^ 2, t43 + t255, t228 + (-qJD(5) + t251) * t86, t98, t21 * t251 - t41 * t86 + t158, t20 * t251 + t41 * t84 - t167, -t254, t250, t249, t247, t98 -(-t147 * t13 - t227) * t97 + (t150 * t98 - t199 * t97 - t86 * t30) * pkin(5) + t246 (-t14 * t97 - t2) * t147 + (t13 * t97 - t185) * t150 + (-t147 * t98 + t169 * t86 - t198 * t97) * pkin(5) + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, t250, t249, t247, t98, t4 * t97 + t246, -t147 * t2 - t150 * t185 + t3 * t97 + t248;];
tauc_reg  = t1;
