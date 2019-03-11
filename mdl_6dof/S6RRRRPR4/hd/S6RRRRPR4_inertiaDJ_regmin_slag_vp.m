% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:20
% EndTime: 2019-03-09 22:09:27
% DurationCPUTime: 2.50s
% Computational Cost: add. (5998->282), mult. (13658->477), div. (0->0), fcn. (13325->10), ass. (0->189)
t191 = sin(qJ(3));
t192 = sin(qJ(2));
t195 = cos(qJ(3));
t196 = cos(qJ(2));
t160 = t191 * t196 + t195 * t192;
t194 = cos(qJ(4));
t183 = qJD(4) * t194;
t226 = t160 * t183;
t159 = t191 * t192 - t195 * t196;
t264 = qJD(2) + qJD(3);
t116 = t264 * t159;
t190 = sin(qJ(4));
t245 = t190 * t116;
t200 = t226 - t245;
t179 = -pkin(2) * t196 - pkin(1);
t106 = t159 * pkin(3) - t160 * pkin(9) + t179;
t263 = -pkin(8) - pkin(7);
t169 = t263 * t192;
t170 = t263 * t196;
t123 = t169 * t191 - t170 * t195;
t119 = t194 * t123;
t242 = t190 * t106 + t119;
t187 = sin(pkin(11));
t188 = cos(pkin(11));
t265 = -t187 * t190 + t188 * t194;
t186 = t194 ^ 2;
t240 = t190 ^ 2 - t186;
t216 = t240 * qJD(4);
t262 = pkin(2) * t195;
t261 = pkin(4) * t187;
t234 = qJD(4) * t190;
t144 = t183 * t188 - t187 * t234;
t260 = pkin(10) * t144;
t153 = t187 * t194 + t188 * t190;
t259 = pkin(10) * t153;
t258 = -qJ(5) - pkin(9);
t117 = t264 * t160;
t204 = qJ(5) * t116 - qJD(5) * t160;
t239 = qJD(2) * t192;
t231 = pkin(2) * t239;
t66 = pkin(3) * t117 + pkin(9) * t116 + t231;
t227 = qJD(2) * t263;
t162 = t192 * t227;
t214 = t196 * t227;
t236 = qJD(3) * t195;
t237 = qJD(3) * t191;
t75 = -t195 * t162 - t169 * t236 - t170 * t237 - t191 * t214;
t223 = t190 * t75 + t194 * t66;
t16 = pkin(4) * t117 + t204 * t194 + (-t119 + (qJ(5) * t160 - t106) * t190) * qJD(4) + t223;
t232 = t106 * t183 + t190 * t66 - t194 * t75;
t20 = -qJ(5) * t226 + (-qJD(4) * t123 + t204) * t190 + t232;
t10 = t187 * t16 + t188 * t20;
t189 = sin(qJ(6));
t193 = cos(qJ(6));
t207 = -t153 * t189 + t193 * t265;
t235 = qJD(4) * t160;
t40 = t116 * t153 - t235 * t265;
t76 = qJD(3) * t123 + t162 * t191 - t195 * t214;
t42 = t200 * pkin(4) + t76;
t25 = -pkin(5) * t40 + t42;
t122 = -t169 * t195 - t191 * t170;
t249 = t160 * t190;
t90 = pkin(4) * t249 + t122;
t98 = t153 * t160;
t58 = pkin(5) * t98 + t90;
t108 = t153 * t193 + t189 * t265;
t143 = t153 * qJD(4);
t62 = qJD(6) * t108 + t193 * t143 + t144 * t189;
t257 = -t207 * t25 + t58 * t62;
t61 = qJD(6) * t207 - t143 * t189 + t144 * t193;
t256 = t25 * t108 + t58 * t61;
t184 = t194 * qJ(5);
t218 = t194 * t106 - t123 * t190;
t48 = pkin(4) * t159 - t160 * t184 + t218;
t57 = -qJ(5) * t249 + t242;
t34 = t187 * t48 + t188 * t57;
t180 = pkin(4) * t234;
t125 = pkin(5) * t143 + t180;
t181 = pkin(2) * t237;
t118 = t125 + t181;
t178 = -pkin(4) * t194 - pkin(3);
t126 = -pkin(5) * t265 + t178;
t124 = t126 - t262;
t255 = -t118 * t207 + t124 * t62;
t254 = t118 * t108 + t124 * t61;
t253 = -t125 * t207 + t126 * t62;
t252 = t125 * t108 + t126 * t61;
t251 = t122 * t183 + t76 * t190;
t250 = t160 * t116;
t248 = t160 * t194;
t244 = t194 * t116;
t176 = pkin(2) * t191 + pkin(9);
t243 = -qJ(5) - t176;
t182 = t194 * qJD(5);
t228 = pkin(2) * t236;
t215 = t194 * t228;
t217 = qJD(4) * t243;
t114 = t190 * t217 + t182 + t215;
t115 = (-qJD(5) - t228) * t190 + t194 * t217;
t68 = t188 * t114 + t187 * t115;
t222 = qJD(4) * t258;
t139 = t190 * t222 + t182;
t140 = -qJD(5) * t190 + t194 * t222;
t93 = t188 * t139 + t187 * t140;
t149 = t243 * t190;
t150 = t176 * t194 + t184;
t103 = t187 * t149 + t188 * t150;
t167 = t258 * t190;
t168 = pkin(9) * t194 + t184;
t121 = t187 * t167 + t188 * t168;
t177 = -pkin(3) - t262;
t241 = t177 * t183 + t190 * t181;
t238 = qJD(2) * t196;
t233 = -0.2e1 * pkin(1) * qJD(2);
t230 = pkin(3) * t234;
t229 = pkin(3) * t183;
t225 = t190 * t183;
t33 = -t187 * t57 + t188 * t48;
t9 = t188 * t16 - t187 * t20;
t224 = t10 * t265 - t34 * t143 - t33 * t144 - t9 * t153;
t102 = t188 * t149 - t150 * t187;
t67 = -t114 * t187 + t188 * t115;
t221 = -t102 * t144 - t103 * t143 - t67 * t153 + t265 * t68;
t120 = t188 * t167 - t168 * t187;
t92 = -t139 * t187 + t188 * t140;
t220 = -t120 * t144 - t121 * t143 - t92 * t153 + t265 * t93;
t219 = -0.4e1 * t190 * t248;
t99 = t265 * t160;
t23 = pkin(5) * t159 - pkin(10) * t99 + t33;
t24 = -pkin(10) * t98 + t34;
t213 = t189 * t24 - t193 * t23;
t212 = t189 * t23 + t193 * t24;
t83 = t102 - t259;
t148 = t265 * pkin(10);
t84 = t148 + t103;
t211 = t189 * t84 - t193 * t83;
t210 = t189 * t83 + t193 * t84;
t94 = t120 - t259;
t95 = t148 + t121;
t209 = t189 * t95 - t193 * t94;
t208 = t189 * t94 + t193 * t95;
t55 = t189 * t99 + t193 * t98;
t56 = -t189 * t98 + t193 * t99;
t206 = t159 * t176 - t160 * t177;
t203 = t177 * t234 - t181 * t194;
t174 = pkin(4) * t188 + pkin(5);
t202 = t174 * t189 + t193 * t261;
t201 = -t174 * t193 + t189 * t261;
t199 = t160 * t234 + t244;
t198 = -t117 * t194 + t159 * t234;
t41 = t153 * t235 - t187 * t245 + t188 * t244;
t5 = pkin(5) * t117 + pkin(10) * t41 + t9;
t6 = pkin(10) * t40 + t10;
t2 = -qJD(6) * t212 - t189 * t6 + t193 * t5;
t1 = qJD(6) * t213 - t189 * t5 - t193 * t6;
t197 = -t116 * t177 - t117 * t176 + (-t159 * t195 + t160 * t191) * qJD(3) * pkin(2);
t172 = 0.2e1 * t225;
t166 = t178 - t262;
t163 = t181 + t180;
t158 = -0.2e1 * t216;
t154 = t160 ^ 2;
t138 = t143 * pkin(10);
t131 = t202 * qJD(6);
t130 = t201 * qJD(6);
t112 = t122 * t234;
t91 = (-t143 * t187 - t144 * t188) * pkin(4);
t87 = 0.2e1 * t159 * t117;
t86 = t117 * t190 + t159 * t183;
t74 = -t138 + t93;
t73 = t92 - t260;
t65 = -t160 * t216 - t190 * t244;
t54 = -t138 + t68;
t53 = t67 - t260;
t52 = qJD(4) * t219 + t240 * t116;
t43 = 0.2e1 * t108 * t61;
t38 = t117 * t207 - t159 * t62;
t37 = t108 * t117 + t159 * t61;
t32 = -t242 * qJD(4) + t223;
t31 = t123 * t234 - t232;
t28 = -qJD(6) * t208 - t189 * t74 + t193 * t73;
t27 = qJD(6) * t209 - t189 * t73 - t193 * t74;
t26 = -0.2e1 * t108 * t62 + 0.2e1 * t207 * t61;
t18 = -qJD(6) * t210 - t189 * t54 + t193 * t53;
t17 = qJD(6) * t211 - t189 * t53 - t193 * t54;
t13 = qJD(6) * t56 - t189 * t41 - t193 * t40;
t12 = -qJD(6) * t55 + t189 * t40 - t193 * t41;
t11 = t108 * t12 + t56 * t61;
t3 = -t108 * t13 + t12 * t207 - t55 * t61 - t56 * t62;
t4 = [0, 0, 0, 0.2e1 * t192 * t238, 0.2e1 * (-t192 ^ 2 + t196 ^ 2) * qJD(2), 0, 0, 0, t192 * t233, t196 * t233, -0.2e1 * t250, 0.2e1 * t116 * t159 - 0.2e1 * t117 * t160, 0, 0, 0, 0.2e1 * t117 * t179 + 0.2e1 * t159 * t231, -0.2e1 * t116 * t179 + 0.2e1 * t160 * t231, -0.2e1 * t154 * t225 - 0.2e1 * t186 * t250, -t116 * t219 + 0.2e1 * t154 * t216, 0.2e1 * t117 * t248 - 0.2e1 * t159 * t199, -0.2e1 * t117 * t249 - 0.2e1 * t159 * t200, t87, 0.2e1 * t117 * t218 + 0.2e1 * t122 * t200 + 0.2e1 * t32 * t159 + 0.2e1 * t76 * t249, -0.2e1 * t242 * t117 - 0.2e1 * t199 * t122 + 0.2e1 * t31 * t159 + 0.2e1 * t76 * t248, -0.2e1 * t10 * t98 + 0.2e1 * t33 * t41 + 0.2e1 * t34 * t40 - 0.2e1 * t9 * t99, 0.2e1 * t10 * t34 + 0.2e1 * t33 * t9 + 0.2e1 * t42 * t90, 0.2e1 * t56 * t12, -0.2e1 * t12 * t55 - 0.2e1 * t13 * t56, 0.2e1 * t117 * t56 + 0.2e1 * t12 * t159, -0.2e1 * t117 * t55 - 0.2e1 * t13 * t159, t87, -0.2e1 * t117 * t213 + 0.2e1 * t58 * t13 + 0.2e1 * t2 * t159 + 0.2e1 * t25 * t55, 0.2e1 * t1 * t159 - 0.2e1 * t117 * t212 + 0.2e1 * t58 * t12 + 0.2e1 * t25 * t56; 0, 0, 0, 0, 0, t238, -t239, 0, -pkin(7) * t238, pkin(7) * t239, 0, 0, -t116, -t117, 0, -t76, t75, t65, t52, t86, -t198, 0, t112 + (-qJD(4) * t206 - t76) * t194 + t197 * t190, t194 * t197 + t206 * t234 + t251, t102 * t41 + t103 * t40 - t67 * t99 - t68 * t98 + t224, t10 * t103 + t102 * t9 + t163 * t90 + t166 * t42 + t33 * t67 + t34 * t68, t11, t3, t37, t38, 0, -t117 * t211 + t118 * t55 + t124 * t13 + t18 * t159 + t257, -t117 * t210 + t118 * t56 + t124 * t12 + t17 * t159 + t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t181, -0.2e1 * t228, t172, t158, 0, 0, 0, 0.2e1 * t203, 0.2e1 * t241, 0.2e1 * t221, 0.2e1 * t102 * t67 + 0.2e1 * t103 * t68 + 0.2e1 * t163 * t166, t43, t26, 0, 0, 0, 0.2e1 * t255, 0.2e1 * t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t117, 0, -t76, t75, t65, t52, t86, -t198, 0, t112 + (pkin(3) * t116 - pkin(9) * t117) * t190 + (-t76 + (-pkin(3) * t160 - pkin(9) * t159) * qJD(4)) * t194, pkin(3) * t199 + pkin(9) * t198 + t251, t120 * t41 + t121 * t40 - t92 * t99 - t93 * t98 + t224, t10 * t121 + t120 * t9 + t178 * t42 + t180 * t90 + t33 * t92 + t34 * t93, t11, t3, t37, t38, 0, -t117 * t209 + t125 * t55 + t126 * t13 + t28 * t159 + t257, -t117 * t208 + t126 * t12 + t125 * t56 + t27 * t159 + t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t228, t172, t158, 0, 0, 0, t203 - t230, -t229 + t241, t220 + t221, t102 * t92 + t103 * t93 + t120 * t67 + t121 * t68 + t163 * t178 + t166 * t180, t43, t26, 0, 0, 0, t253 + t255, t252 + t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t158, 0, 0, 0, -0.2e1 * t230, -0.2e1 * t229, 0.2e1 * t220, 0.2e1 * t120 * t92 + 0.2e1 * t121 * t93 + 0.2e1 * t178 * t180, t43, t26, 0, 0, 0, 0.2e1 * t253, 0.2e1 * t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, -t200, t117, t32, t31 (t187 * t40 + t188 * t41) * pkin(4) (t10 * t187 + t188 * t9) * pkin(4), 0, 0, t12, -t13, t117, -t117 * t201 - t131 * t159 + t2, -t117 * t202 + t130 * t159 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, -t234, 0, -t176 * t183 - t190 * t228, t176 * t234 - t215, t91 (t187 * t68 + t188 * t67) * pkin(4), 0, 0, t61, -t62, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, -t234, 0, -pkin(9) * t183, pkin(9) * t234, t91 (t187 * t93 + t188 * t92) * pkin(4), 0, 0, t61, -t62, 0, t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t131, 0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, 0, 0, 0, 0, 0, t62, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, 0, 0, 0, 0, 0, t62, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, t117, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62, 0, t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
