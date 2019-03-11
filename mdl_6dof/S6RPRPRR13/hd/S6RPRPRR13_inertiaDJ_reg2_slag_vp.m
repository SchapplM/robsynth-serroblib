% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR13_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:21
% EndTime: 2019-03-09 04:25:41
% DurationCPUTime: 7.81s
% Computational Cost: add. (10998->410), mult. (31795->737), div. (0->0), fcn. (33809->12), ass. (0->212)
t250 = sin(qJ(3));
t89 = sin(pkin(7));
t188 = t89 * t250;
t223 = cos(pkin(7));
t91 = cos(pkin(12));
t176 = t91 * t223;
t224 = cos(pkin(6));
t178 = t89 * t224;
t251 = cos(qJ(3));
t88 = sin(pkin(12));
t90 = sin(pkin(6));
t53 = t250 * t178 + (t250 * t176 + t251 * t88) * t90;
t236 = t90 * t91;
t60 = -t224 * t223 + t89 * t236;
t271 = (t60 * t188 + t223 * t53) * qJD(3);
t111 = qJD(3) * t53;
t270 = -0.2e1 * t111 * t60;
t187 = t90 * t250;
t177 = t90 * t223;
t150 = t251 * t177;
t232 = t91 * t150 + t251 * t178;
t136 = t88 * t187 - t232;
t269 = 0.2e1 * t136 * t111;
t183 = pkin(1) * t224;
t129 = t224 * pkin(2) + t91 * t183;
t117 = (-t223 * pkin(9) - qJ(2)) * t90 * t88 + t129;
t190 = -pkin(2) * t91 - pkin(1);
t134 = (-pkin(9) * t88 * t89 + t190) * t90;
t268 = t117 * t223 + t89 * t134;
t189 = t89 * t251;
t49 = t136 * qJD(3);
t180 = qJD(3) * t250;
t80 = t89 * t180;
t266 = -t111 * t188 + t49 * t189 + t53 * t80;
t265 = -0.2e1 * t111 * t53 + 0.2e1 * t49 * t136;
t95 = cos(qJ(5));
t186 = t95 * t251;
t93 = sin(qJ(5));
t130 = -t89 * t186 - t93 * t223;
t210 = t130 * qJD(5);
t96 = -pkin(3) - pkin(10);
t163 = pkin(5) * t95 + pkin(11) * t93;
t94 = cos(qJ(6));
t260 = t94 * t163;
t114 = t88 * t183 + qJ(2) * t236 + (t90 * t176 + t178) * pkin(9);
t30 = t251 * t114 + t250 * t268;
t25 = t60 * qJ(4) - t30;
t20 = -t136 * pkin(4) - t25;
t132 = t136 * t95;
t41 = -t93 * t60 - t132;
t42 = t136 * t93 - t60 * t95;
t101 = t41 * pkin(5) - t42 * pkin(11) + t20;
t29 = -t250 * t114 + t251 * t268;
t26 = t60 * pkin(3) - t29;
t103 = t53 * pkin(4) + t60 * pkin(10) + t26;
t222 = qJD(2) * t90;
t203 = t88 * t222;
t170 = t89 * t203;
t27 = pkin(3) * t111 + t49 * qJ(4) - t53 * qJD(4) + t170;
t257 = -pkin(10) * t111 - qJD(5) * t103 - t27;
t226 = t53 * qJ(4);
t106 = t89 * t117 - t223 * t134 + t96 * t136 + t226;
t23 = (t150 * t88 + t187 * t91) * qJD(2) + t30 * qJD(3);
t258 = -t49 * pkin(4) + qJD(5) * t106 + t23;
t5 = t257 * t95 - t258 * t93;
t259 = t49 * pkin(11) - qJD(6) * t101 + t5;
t164 = t93 * t180;
t116 = t89 * t164 + t210;
t216 = qJD(5) * t95;
t167 = t93 * t189;
t63 = t95 * t223 - t167;
t256 = t116 * t93 + t63 * t216;
t255 = -qJD(6) * t188 - t116;
t92 = sin(qJ(6));
t84 = t92 ^ 2;
t86 = t94 ^ 2;
t230 = t84 - t86;
t175 = qJD(6) * t230;
t253 = 0.2e1 * qJD(4);
t252 = pkin(5) * t93;
t218 = qJD(5) * t93;
t31 = -qJD(5) * t132 - t111 * t93 - t60 * t218;
t34 = t42 * t94 + t53 * t92;
t15 = t34 * qJD(6) - t31 * t92 + t49 * t94;
t249 = t15 * t92;
t248 = t15 * t94;
t214 = qJD(6) * t92;
t16 = -t49 * t92 - t42 * t214 + (qJD(6) * t53 - t31) * t94;
t247 = t16 * t92;
t246 = t23 * t53;
t245 = t23 * t60;
t244 = t31 * t95;
t32 = qJD(5) * t42 - t111 * t95;
t243 = t32 * t93;
t33 = t42 * t92 - t53 * t94;
t242 = t33 * t92;
t241 = t33 * t94;
t240 = t34 * t92;
t239 = t34 * t94;
t238 = t60 * t49;
t173 = qJD(5) * t223;
t54 = -qJD(5) * t167 + (t173 - t80) * t95;
t237 = t130 * t54;
t235 = t93 * t96;
t234 = t94 * t95;
t233 = t95 * t96;
t13 = t93 * t103 - t95 * t106;
t181 = qJD(3) * t251;
t165 = t89 * t181;
t231 = qJ(4) * t165 + qJD(4) * t188;
t229 = t84 + t86;
t85 = t93 ^ 2;
t87 = t95 ^ 2;
t228 = t85 - t87;
t227 = t85 + t87;
t225 = t89 * qJ(2);
t221 = qJD(4) * t60;
t220 = qJD(5) * t33;
t219 = qJD(5) * t34;
t217 = qJD(5) * t94;
t215 = qJD(5) * t96;
t213 = qJD(6) * t94;
t212 = qJD(6) * t95;
t211 = qJD(6) * t96;
t209 = qJ(4) * qJD(5);
t208 = 0.2e1 * t41 * t32;
t38 = -0.2e1 * t53 * t49;
t207 = -0.2e1 * pkin(5) * qJD(6);
t206 = t92 * t235;
t205 = t92 * t233;
t204 = t94 * t235;
t202 = qJD(5) * t41 * t92;
t201 = t41 * t217;
t200 = t93 * t217;
t199 = t93 * t215;
t198 = t95 * t215;
t197 = t92 * t212;
t196 = t92 * t211;
t195 = t94 * t212;
t194 = t93 * t210;
t192 = t92 * t213;
t191 = t93 * t216;
t185 = t251 * t23;
t182 = t95 * t229;
t179 = qJD(5) * t250;
t174 = t228 * qJD(5);
t78 = 0.2e1 * t191;
t172 = t92 * t200;
t171 = t87 * t192;
t162 = -pkin(11) * t95 + t252;
t11 = pkin(11) * t53 + t13;
t7 = t101 * t94 - t92 * t11;
t8 = t101 * t92 + t94 * t11;
t161 = t7 * t94 + t8 * t92;
t160 = t7 * t92 - t8 * t94;
t159 = t240 + t241;
t55 = t94 * t188 - t92 * t63;
t56 = t92 * t188 + t94 * t63;
t158 = t55 * t94 + t56 * t92;
t157 = t55 * t92 - t56 * t94;
t151 = qJ(4) + t162;
t135 = t94 * t151;
t58 = t135 - t206;
t59 = t92 * t151 + t204;
t156 = t58 * t94 + t59 * t92;
t155 = t58 * t92 - t59 * t94;
t153 = -0.2e1 * t224 * t222;
t152 = qJD(2) * t88 * t177;
t12 = t103 * t95 + t106 * t93;
t10 = -t53 * pkin(5) - t12;
t6 = t257 * t93 + t258 * t95;
t4 = t49 * pkin(5) - t6;
t147 = t10 * t213 + t4 * t92;
t146 = t10 * t214 - t4 * t94;
t145 = t41 * t216 + t243;
t144 = -t53 * t216 + t49 * t93;
t143 = -t54 * t95 - t194;
t142 = t41 * t213 + t32 * t92;
t141 = t41 * t214 - t32 * t94;
t140 = -t130 * t213 + t54 * t92;
t139 = -t130 * t214 - t54 * t94;
t138 = t223 * t190;
t137 = 0.2e1 * (t88 ^ 2 + t91 ^ 2) * t90 ^ 2 * qJD(2);
t65 = -t197 - t200;
t67 = t92 * t218 - t195;
t133 = t89 * t136;
t128 = t60 * t165 - t223 * t49;
t126 = t89 * t129;
t22 = -t251 * t91 * t222 - t29 * qJD(3) + t250 * t152;
t21 = t22 + t221;
t17 = -pkin(4) * t111 - t21;
t98 = t32 * pkin(5) + t31 * pkin(11) + t17;
t1 = t11 * t214 + t259 * t94 - t92 * t98;
t2 = -t11 * t213 + t259 * t92 + t94 * t98;
t122 = -t161 * qJD(6) - t1 * t94 - t2 * t92;
t121 = -t5 * t93 + t6 * t95 + (-t12 * t93 + t13 * t95) * qJD(5);
t120 = -t248 + t247 + (t239 + t242) * qJD(6);
t36 = -t92 * t165 + t63 * t214 + t255 * t94;
t37 = t94 * t165 - t63 * t213 + t255 * t92;
t119 = -t158 * qJD(6) - t36 * t94 - t37 * t92;
t43 = t93 * t196 - t92 * (t163 * qJD(5) + qJD(4)) - qJD(6) * t135 - t94 * t198;
t44 = t94 * qJD(4) - t59 * qJD(6) + (-t205 + t260) * qJD(5);
t118 = -t156 * qJD(6) - t43 * t94 - t44 * t92;
t82 = qJ(4) * t253;
t66 = t93 * t213 + t92 * t216;
t64 = t93 * t214 - t94 * t216;
t57 = t95 * t175 + t172;
t39 = -t126 + (t88 * t225 + t138) * t90;
t35 = -t53 * t218 - t49 * t95;
t28 = t143 + t256;
t24 = -t126 - t226 - t232 * pkin(3) + (t138 + (t250 * pkin(3) + t225) * t88) * t90;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t153, t91 * t153, t137, qJ(2) * t137, t38, t265, 0.2e1 * t238, t269, -t270, 0, 0.2e1 * t111 * t39 + 0.2e1 * t133 * t203 + 0.2e1 * t245, 0.2e1 * t53 * t170 - 0.2e1 * t22 * t60 - 0.2e1 * t39 * t49, -0.2e1 * t111 * t30 + 0.2e1 * t22 * t136 + 0.2e1 * t29 * t49 + 0.2e1 * t246, 0.2e1 * t39 * t170 - 0.2e1 * t22 * t30 - 0.2e1 * t23 * t29, 0, -0.2e1 * t238, t270, t38, t265, t269, 0.2e1 * t111 * t25 + 0.2e1 * t136 * t21 - 0.2e1 * t26 * t49 + 0.2e1 * t246, -0.2e1 * t111 * t24 - 0.2e1 * t136 * t27 - 0.2e1 * t245, 0.2e1 * t21 * t60 + 0.2e1 * t24 * t49 - 0.2e1 * t27 * t53, 0.2e1 * t21 * t25 + 0.2e1 * t23 * t26 + 0.2e1 * t24 * t27, -0.2e1 * t42 * t31, 0.2e1 * t31 * t41 - 0.2e1 * t32 * t42, -0.2e1 * t31 * t53 - 0.2e1 * t42 * t49, t208, -0.2e1 * t32 * t53 + 0.2e1 * t41 * t49, t38, -0.2e1 * t12 * t49 + 0.2e1 * t17 * t41 + 0.2e1 * t20 * t32 + 0.2e1 * t53 * t6, 0.2e1 * t13 * t49 + 0.2e1 * t17 * t42 - 0.2e1 * t20 * t31 + 0.2e1 * t5 * t53, 0.2e1 * t12 * t31 - 0.2e1 * t13 * t32 + 0.2e1 * t41 * t5 - 0.2e1 * t42 * t6, 0.2e1 * t12 * t6 - 0.2e1 * t13 * t5 + 0.2e1 * t17 * t20, 0.2e1 * t34 * t16, -0.2e1 * t15 * t34 - 0.2e1 * t16 * t33, 0.2e1 * t16 * t41 + 0.2e1 * t32 * t34, 0.2e1 * t33 * t15, -0.2e1 * t15 * t41 - 0.2e1 * t32 * t33, t208, 0.2e1 * t10 * t15 + 0.2e1 * t2 * t41 + 0.2e1 * t32 * t7 + 0.2e1 * t33 * t4, 0.2e1 * t1 * t41 + 0.2e1 * t10 * t16 - 0.2e1 * t32 * t8 + 0.2e1 * t34 * t4, 0.2e1 * t1 * t33 - 0.2e1 * t15 * t8 - 0.2e1 * t16 * t7 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t8 + 0.2e1 * t10 * t4 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, t128, -t136 * t165 + t266 (t152 - t250 * t22 - t185 + (-t250 * t29 + t251 * t30) * qJD(3)) * t89, 0, 0, 0, 0, 0, 0, -t133 * t181 + t266, -t271, -t128, t27 * t223 + (-t250 * t21 - t185 + (-t251 * t25 + t250 * t26) * qJD(3)) * t89, 0, 0, 0, 0, 0, 0, -t130 * t49 - t54 * t53 + (t41 * t181 + t250 * t32) * t89, -t116 * t53 + t165 * t42 - t188 * t31 + t63 * t49, -t116 * t41 + t130 * t31 - t63 * t32 + t54 * t42, t116 * t13 - t12 * t54 + t130 * t6 + t165 * t20 + t17 * t188 - t5 * t63, 0, 0, 0, 0, 0, 0, -t130 * t15 + t32 * t55 + t33 * t54 + t37 * t41, -t130 * t16 - t32 * t56 + t34 * t54 + t36 * t41, -t15 * t56 - t16 * t55 + t33 * t36 - t34 * t37, -t1 * t56 + t10 * t54 - t130 * t4 + t2 * t55 - t36 * t8 + t37 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63 * t93 * t173 - 0.2e1 * t237 + 0.2e1 * (t63 * (-qJD(5) * t186 + t164) + t251 * t80) * t89, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t36 * t56 + 0.2e1 * t37 * t55 - 0.2e1 * t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t111, 0, -t23, t22, 0, 0, 0, t49, t111, 0, 0, 0, t49 * pkin(3) - qJ(4) * t111 - qJD(4) * t136, t23, -t22 - 0.2e1 * t221, -pkin(3) * t23 - qJ(4) * t21 - qJD(4) * t25, -t218 * t42 - t244, t31 * t93 - t32 * t95 + (t41 * t93 - t42 * t95) * qJD(5), t35, t145, t144, 0, -t49 * t233 + qJ(4) * t32 + qJD(4) * t41 + t17 * t93 + (t20 * t95 - t53 * t235) * qJD(5), t49 * t235 - qJ(4) * t31 + qJD(4) * t42 + t17 * t95 + (-t20 * t93 - t53 * t233) * qJD(5) (t31 * t96 - t6) * t95 + (-t32 * t96 + t5) * t93 + ((-t41 * t96 - t13) * t95 + (t42 * t96 + t12) * t93) * qJD(5), t17 * qJ(4) + t20 * qJD(4) + t121 * t96, t16 * t234 + t34 * t65, t159 * t218 + (-t248 - t247 + (-t239 + t242) * qJD(6)) * t95 (t16 - t201) * t93 + (-t141 + t219) * t95, t95 * t249 - t33 * t67 (-t15 + t202) * t93 + (-t142 - t220) * t95, t145, t58 * t32 + t44 * t41 + (t2 + (-t10 * t92 + t33 * t96) * qJD(5)) * t93 + (qJD(5) * t7 - t15 * t96 + t147) * t95, -t59 * t32 + t43 * t41 + (t1 + (-t10 * t94 + t34 * t96) * qJD(5)) * t93 + (-qJD(5) * t8 - t16 * t96 - t146) * t95, -t15 * t59 - t16 * t58 + t33 * t43 - t34 * t44 + t161 * t218 + (qJD(6) * t160 + t1 * t92 - t2 * t94) * t95, -t1 * t59 + t2 * t58 - t8 * t43 + t7 * t44 + (t10 * t218 - t4 * t95) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t165, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t165, -pkin(3) * t80 + t231, 0, 0, 0, 0, 0, 0 (t179 * t95 + t181 * t93) * t89 (-t179 * t93 + t181 * t95) * t89, -t28, -t54 * t233 + t231 + (-t194 + t256) * t96, 0, 0, 0, 0, 0, 0 (t210 * t92 + t37) * t93 + (qJD(5) * t55 + t140) * t95 (t210 * t94 + t36) * t93 + (-qJD(5) * t56 - t139) * t95, t158 * t218 + (qJD(6) * t157 + t36 * t92 - t37 * t94) * t95, t143 * t96 - t36 * t59 + t37 * t58 - t56 * t43 + t55 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t82, -0.2e1 * t191, 0.2e1 * t174, 0, t78, 0, 0, 0.2e1 * qJD(4) * t93 + 0.2e1 * t209 * t95, 0.2e1 * qJD(4) * t95 - 0.2e1 * t209 * t93, 0, t82, -0.2e1 * t191 * t86 - 0.2e1 * t171, 0.4e1 * t172 * t95 + 0.2e1 * t175 * t87, -0.2e1 * t93 * t197 - 0.2e1 * t228 * t217, -0.2e1 * t191 * t84 + 0.2e1 * t171, 0.2e1 * t174 * t92 - 0.2e1 * t195 * t93, t78, -0.2e1 * t87 * t94 * t211 + 0.2e1 * t44 * t93 + 0.2e1 * (t58 + 0.2e1 * t206) * t216, 0.2e1 * t87 * t196 + 0.2e1 * t43 * t93 + 0.2e1 * (-t59 + 0.2e1 * t204) * t216, 0.2e1 * t156 * t218 + 0.2e1 * (qJD(6) * t155 + t43 * t92 - t44 * t94) * t95, -0.2e1 * t191 * t96 ^ 2 - 0.2e1 * t59 * t43 + 0.2e1 * t58 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, t23, 0, 0, 0, 0, 0, 0, t35, t144, t244 - t243 + (-t41 * t95 + t42 * t93) * qJD(5), t121, 0, 0, 0, 0, 0, 0 (-t15 - t202) * t95 + (-t142 + t220) * t93 (-t16 - t201) * t95 + (t141 + t219) * t93 (t240 - t241) * t216 + t120 * t93 (-qJD(5) * t160 - t4) * t95 + (qJD(5) * t10 + t122) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-qJD(5) * t157 - t54) * t95 + (t119 - t210) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227 * t213, t227 * t214, 0, -t155 * t216 + (t118 - 0.2e1 * t198) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t229) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t32, -t49, t6, t5, 0, 0, t213 * t34 + t247, -qJD(6) * t159 + t16 * t94 - t249, t142, t214 * t33 - t248, -t141, 0, -pkin(5) * t15 - pkin(11) * t142 + t146, -pkin(5) * t16 + pkin(11) * t141 + t147, pkin(11) * t120 + t122, -pkin(5) * t4 + pkin(11) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t116, 0, 0, 0, 0, 0, 0, 0, 0, t139, t140, t119, -pkin(5) * t54 + pkin(11) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, 0, -t216, 0, -t199, -t198, 0, 0, -t57, -0.4e1 * t95 * t192 + t230 * t218, t66, t57, -t64, 0 (-t205 - t260) * qJD(6) + (t162 * t92 - t204) * qJD(5) (t163 * t92 - t94 * t233) * qJD(6) + (-pkin(11) * t234 + (pkin(5) * t94 + t92 * t96) * t93) * qJD(5), t118, -pkin(5) * t199 + pkin(11) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, -t216, 0, 0, 0, 0, 0, 0, 0, 0, t65, t67, qJD(5) * t182 (pkin(11) * t182 - t252) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t192, -0.2e1 * t175, 0, -0.2e1 * t192, 0, 0, t92 * t207, t94 * t207, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t15, t32, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t67, t216, t44, t43, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, 0, -t214, 0, -pkin(11) * t213, pkin(11) * t214, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
