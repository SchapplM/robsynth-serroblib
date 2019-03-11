% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:41
% EndTime: 2019-03-09 04:44:51
% DurationCPUTime: 3.26s
% Computational Cost: add. (4935->375), mult. (12648->444), div. (0->0), fcn. (9252->6), ass. (0->181)
t144 = sin(pkin(9));
t145 = cos(pkin(9));
t147 = sin(qJ(3));
t149 = cos(qJ(3));
t265 = -t144 * t147 + t149 * t145;
t268 = t265 * qJD(1);
t109 = qJD(4) - t268;
t148 = cos(qJ(4));
t146 = sin(qJ(4));
t201 = qJD(3) * t146;
t126 = t144 * t149 + t145 * t147;
t251 = t126 * qJD(1);
t89 = t148 * t251 + t201;
t221 = t89 * t109;
t199 = qJD(4) * t148;
t49 = t251 * t199 + (qJD(4) + t268) * t201;
t269 = t49 - t221;
t120 = t126 * qJD(3);
t103 = qJD(1) * t120;
t267 = t103 * qJ(5) + t109 * qJD(5);
t239 = pkin(7) + qJ(2);
t131 = t239 * t144;
t127 = qJD(1) * t131;
t132 = t239 * t145;
t128 = qJD(1) * t132;
t77 = -t147 * t127 + t149 * t128;
t266 = qJD(3) * t77;
t264 = t251 * qJD(3);
t196 = t148 * qJD(3);
t87 = t146 * t251 - t196;
t262 = t49 * qJ(6) + t87 * qJD(6);
t102 = t109 * qJ(5);
t139 = -pkin(2) * t145 - pkin(1);
t129 = t139 * qJD(1) + qJD(2);
t54 = -pkin(3) * t268 - pkin(8) * t251 + t129;
t72 = qJD(3) * pkin(8) + t77;
t29 = t146 * t54 + t148 * t72;
t22 = t102 + t29;
t232 = t109 * t22;
t200 = qJD(4) * t146;
t162 = t265 * qJD(2);
t254 = -t149 * t127 - t147 * t128;
t39 = qJD(1) * t162 + qJD(3) * t254;
t119 = t265 * qJD(3);
t158 = qJD(1) * t119;
t59 = t103 * pkin(3) - pkin(8) * t158;
t187 = t146 * t39 - t148 * t59 + t72 * t199 + t54 * t200;
t99 = t103 * pkin(4);
t7 = -t99 + t187;
t261 = -t7 + t232;
t184 = qJD(3) * pkin(3) + t254;
t160 = t89 * qJ(5) + t184;
t246 = pkin(4) + pkin(5);
t17 = -t246 * t87 + qJD(6) + t160;
t260 = (qJD(6) + t17) * t89;
t259 = 0.2e1 * t267;
t104 = t109 ^ 2;
t86 = t89 ^ 2;
t258 = -t104 - t86;
t257 = -t146 * qJD(5) - t77;
t230 = t109 * t87;
t48 = -qJD(4) * t196 - t148 * t158 + t200 * t251;
t157 = t48 - t230;
t83 = -t131 * t147 + t132 * t149;
t256 = t83 * qJD(3);
t28 = -t146 * t72 + t148 * t54;
t204 = qJD(5) - t28;
t255 = t120 * qJ(5) - qJD(5) * t265;
t214 = t268 * t146;
t96 = t148 * t103;
t172 = t96 + (-t200 + t214) * t109;
t229 = t251 * t87;
t154 = t172 - t229;
t253 = -t149 * t131 - t147 * t132;
t18 = t89 * qJ(6) + t28;
t205 = qJD(5) - t18;
t12 = -t246 * t109 + t205;
t169 = -t146 * t59 - t148 * t39 - t54 * t199 + t72 * t200;
t6 = -t169 + t267;
t2 = t6 + t262;
t250 = t109 * t12 + t2;
t245 = pkin(8) * t103;
t30 = pkin(4) * t87 - t160;
t249 = t109 * t30 - t245;
t208 = t146 * qJ(5);
t248 = -t246 * t148 - t208;
t247 = t87 ^ 2;
t195 = qJD(1) * qJD(2);
t40 = t126 * t195 + t266;
t9 = t49 * pkin(4) + t48 * qJ(5) - t89 * qJD(5) + t40;
t5 = -pkin(5) * t49 - t9;
t244 = t5 * t146;
t243 = t5 * t148;
t242 = t89 * t87;
t241 = t9 * t146;
t240 = t9 * t148;
t238 = pkin(8) - qJ(6);
t216 = qJ(5) * t148;
t170 = -t246 * t146 + t216;
t237 = -t109 * t170 + t257;
t236 = -t146 * t49 - t87 * t199;
t73 = pkin(3) * t251 - pkin(8) * t268;
t235 = t146 * t73 + t148 * t254;
t75 = -pkin(3) * t265 - pkin(8) * t126 + t139;
t234 = t146 * t75 + t148 * t83;
t233 = qJ(5) * t49;
t231 = t109 * t29;
t19 = qJ(6) * t87 + t29;
t14 = t102 + t19;
t228 = t14 * t109;
t227 = t146 * t89;
t226 = t148 * t87;
t225 = t40 * t146;
t224 = t40 * t148;
t223 = t48 * t146;
t222 = t87 * qJ(5);
t220 = t89 * t251;
t24 = qJ(5) * t251 + t235;
t219 = -qJ(6) * t214 - t148 * qJD(6) - t238 * t200 - t24;
t180 = pkin(4) * t146 - t216;
t218 = t109 * t180 + t257;
t134 = t238 * t148;
t69 = t146 * t254;
t217 = qJD(4) * t134 - t146 * qJD(6) - t69 - (-qJ(6) * t268 - t73) * t148 + t246 * t251;
t215 = qJ(6) * t126;
t185 = t109 * t148;
t213 = t119 * t146;
t212 = t119 * t148;
t94 = t146 * t103;
t202 = t144 ^ 2 + t145 ^ 2;
t198 = qJD(5) * t148;
t55 = t253 * qJD(3) + t162;
t194 = t146 * t55 + t83 * t199 + t75 * t200;
t74 = pkin(3) * t120 - pkin(8) * t119;
t193 = t146 * t74 + t148 * t55 + t75 * t199;
t31 = -qJ(5) * t265 + t234;
t192 = pkin(8) * t200;
t191 = t126 * t200;
t190 = t126 * t199;
t79 = t146 * t83;
t188 = t148 * t75 - t79;
t186 = t202 * qJD(1) ^ 2;
t165 = t48 * qJ(6) + t7;
t155 = -t103 * pkin(5) + t165;
t1 = -t89 * qJD(6) + t155;
t183 = -t1 + t228;
t21 = -pkin(4) * t109 + t204;
t182 = t109 * t21 + t6;
t181 = pkin(4) * t148 + t208;
t179 = -t146 * t22 + t148 * t21;
t178 = t226 + t227;
t177 = t148 * t74 - t194;
t175 = -qJ(6) * t119 - qJD(6) * t126;
t174 = t109 * t199 - t185 * t268 + t94;
t173 = 0.2e1 * t202 * t195;
t171 = -t148 * t48 - t89 * t200;
t168 = -t83 * t200 + t193;
t167 = -t109 * t184 - t245;
t166 = t30 * t89 + t7;
t159 = t174 + t220;
t156 = t109 * t28 + t169;
t153 = t242 - t264;
t56 = t126 * qJD(2) + t256;
t133 = t238 * t146;
t130 = -pkin(3) - t181;
t123 = pkin(3) - t248;
t41 = pkin(4) * t89 + t222;
t36 = t180 * t126 - t253;
t34 = -t246 * t89 - t222;
t33 = t170 * t126 + t253;
t32 = pkin(4) * t265 - t188;
t26 = -pkin(4) * t251 - t148 * t73 + t69;
t23 = t146 * t215 + t31;
t16 = t79 + (-t75 - t215) * t148 + t246 * t265;
t13 = t180 * t119 + (t181 * qJD(4) - t198) * t126 + t56;
t11 = -t120 * pkin(4) - t177;
t10 = -t256 + t170 * t119 + (t248 * qJD(4) - qJD(2) + t198) * t126;
t8 = t168 + t255;
t4 = qJ(6) * t190 + (-qJD(4) * t83 - t175) * t146 + t193 + t255;
t3 = qJ(6) * t191 - t246 * t120 + (t175 - t74) * t148 + t194;
t15 = [0, 0, 0, 0, 0, t173, qJ(2) * t173, t119 * t251 + t126 * t158, -t126 * t103 + t119 * t268 - t120 * t251 + t158 * t265, t119 * qJD(3), -t120 * qJD(3), 0, -qJD(3) * t56 + t103 * t139 + t120 * t129, t129 * t119 + (t139 * t268 - t55) * qJD(3), t171 * t126 + t89 * t212, -t178 * t119 + (t223 - t148 * t49 + (t146 * t87 - t148 * t89) * qJD(4)) * t126, t126 * t96 + t89 * t120 + t48 * t265 + (-t191 + t212) * t109, -t126 * t94 - t87 * t120 + t49 * t265 + (-t190 - t213) * t109, -t103 * t265 + t109 * t120, t177 * t109 + t188 * t103 + t187 * t265 + t28 * t120 + t56 * t87 - t253 * t49 - t184 * t213 + (-t184 * t199 + t225) * t126, -t168 * t109 - t234 * t103 - t169 * t265 - t29 * t120 + t56 * t89 + t253 * t48 - t184 * t212 + (t184 * t200 + t224) * t126, t30 * t213 - t32 * t103 - t11 * t109 - t21 * t120 + t7 * t265 + t13 * t87 + t36 * t49 + (t30 * t199 + t241) * t126, t11 * t89 - t31 * t49 - t32 * t48 - t8 * t87 + t179 * t119 + (-t6 * t146 + t7 * t148 + (-t146 * t21 - t148 * t22) * qJD(4)) * t126, -t30 * t212 + t31 * t103 + t8 * t109 + t22 * t120 - t6 * t265 - t13 * t89 + t36 * t48 + (t200 * t30 - t240) * t126, t11 * t21 + t13 * t30 + t22 * t8 + t31 * t6 + t32 * t7 + t36 * t9, -t17 * t213 + t1 * t265 - t10 * t87 - t16 * t103 - t3 * t109 - t12 * t120 - t33 * t49 + (-t17 * t199 - t244) * t126, t17 * t212 + t10 * t89 + t23 * t103 + t4 * t109 + t14 * t120 - t2 * t265 - t33 * t48 + (-t17 * t200 + t243) * t126, t16 * t48 + t23 * t49 - t3 * t89 + t4 * t87 + (-t12 * t148 + t14 * t146) * t119 + (-t1 * t148 + t2 * t146 + (t12 * t146 + t14 * t148) * qJD(4)) * t126, t1 * t16 + t10 * t17 + t12 * t3 + t14 * t4 + t2 * t23 + t33 * t5; 0, 0, 0, 0, 0, -t186, -qJ(2) * t186, 0, 0, 0, 0, 0, 0.2e1 * t264, 0.2e1 * t268 * qJD(3), 0, 0, 0, 0, 0, t154, -t104 * t148 - t220 - t94, t154 (t226 - t227) * t268 - t171 + t236, t159, t182 * t146 + t261 * t148 - t251 * t30, t154, t159, t269 * t146 - t157 * t148, t250 * t146 + t183 * t148 + t17 * t251; 0, 0, 0, 0, 0, 0, 0, -t251 * t268, t251 ^ 2 - t268 ^ 2, 0, 0, 0, -t129 * t251 + t266 - t40 (-qJD(2) - t129) * t268, t89 * t185 - t223, t178 * t268 + t171 + t236, t174 - t220, t172 + t229, -t109 * t251, -pkin(3) * t49 - t28 * t251 - t224 - t77 * t87 + (t69 + (-pkin(8) * qJD(4) - t73) * t148) * t109 + t167 * t146, pkin(3) * t48 + t29 * t251 + t225 - t77 * t89 + (t192 + t235) * t109 + t167 * t148, t21 * t251 + t130 * t49 - t240 + t218 * t87 + (-pkin(8) * t199 + t26) * t109 + t249 * t146, t24 * t87 - t26 * t89 + ((qJD(4) * t89 - t49) * pkin(8) + t182) * t148 + ((qJD(4) * t87 - t48) * pkin(8) - t261) * t146, -t22 * t251 + t130 * t48 - t241 - t218 * t89 + (-t24 - t192) * t109 - t249 * t148, t9 * t130 - t21 * t26 - t22 * t24 + t218 * t30 + (qJD(4) * t179 + t7 * t146 + t6 * t148) * pkin(8), -t133 * t103 + t12 * t251 - t123 * t49 + t237 * t87 + t243 + (-t146 * t17 - t217) * t109, t134 * t103 + t219 * t109 - t123 * t48 - t14 * t251 + t17 * t185 - t237 * t89 + t244, t133 * t48 + t134 * t49 + t183 * t146 - t250 * t148 - t217 * t89 + t219 * t87, t1 * t133 + t217 * t12 + t5 * t123 + t2 * t134 + t219 * t14 - t237 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t86 - t247, -t157, -t269, t103, t184 * t89 - t187 + t231, -t184 * t87 + t156, -t41 * t87 - t166 + t231 + t99, pkin(4) * t48 - t233 + (t22 - t29) * t89 + (t21 - t204) * t87, -t30 * t87 + t41 * t89 - t156 + t259, -t7 * pkin(4) + t6 * qJ(5) + t204 * t22 - t21 * t29 - t30 * t41, t19 * t109 + t34 * t87 + t260 + (pkin(5) + t246) * t103 - t165, -t109 * t18 + t17 * t87 - t34 * t89 - t169 + t259 + t262, t233 - t246 * t48 + (-t14 + t19) * t89 + (-t12 + t205) * t87, t2 * qJ(5) - t1 * t246 - t12 * t19 + t14 * t205 - t17 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t157, t258, t166 - t232, t153, t258, t157, t155 - t228 - t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49 - t221, -t48 - t230, -t86 - t247, t12 * t89 - t14 * t87 + t5;];
tauc_reg  = t15;
