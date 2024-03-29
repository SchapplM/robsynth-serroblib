% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:57
% EndTime: 2022-01-20 11:30:59
% DurationCPUTime: 2.30s
% Computational Cost: add. (2646->256), mult. (6000->354), div. (0->0), fcn. (5044->8), ass. (0->207)
t159 = sin(qJ(5));
t155 = t159 ^ 2;
t162 = cos(qJ(5));
t156 = t162 ^ 2;
t146 = t156 - t155;
t203 = -qJD(2) - qJD(3);
t192 = qJD(1) - t203;
t273 = t192 * t146;
t215 = t155 + t156;
t193 = t156 / 0.2e1 + t155 / 0.2e1;
t163 = cos(qJ(3));
t164 = cos(qJ(2));
t253 = t164 * pkin(1);
t195 = pkin(2) + t253;
t139 = t163 * t195;
t160 = sin(qJ(3));
t161 = sin(qJ(2));
t229 = t160 * t161;
t125 = pkin(1) * t229 - t139;
t120 = pkin(3) - t125;
t158 = cos(pkin(9));
t183 = t160 * t195;
t227 = t163 * t161;
t126 = pkin(1) * t227 + t183;
t157 = sin(pkin(9));
t236 = t157 * t126;
t69 = t158 * t120 - t236;
t65 = -pkin(4) - t69;
t272 = -t65 / 0.2e1;
t232 = t158 * t126;
t237 = t157 * t125;
t78 = t232 - t237;
t271 = t78 / 0.2e1;
t79 = -t158 * t125 - t236;
t270 = -t79 / 0.2e1;
t269 = t79 / 0.2e1;
t228 = t160 * t164;
t131 = (t227 + t228) * pkin(1);
t119 = t158 * t131;
t226 = t163 * t164;
t132 = (t226 - t229) * pkin(1);
t234 = t157 * t132;
t84 = t119 + t234;
t268 = -t84 / 0.2e1;
t231 = t158 * t132;
t235 = t157 * t131;
t85 = t231 - t235;
t267 = -t85 / 0.2e1;
t266 = -t232 / 0.2e1;
t90 = t160 * pkin(2);
t148 = t157 * t90;
t254 = t163 * pkin(2);
t151 = pkin(3) + t254;
t121 = t158 * t151 - t148;
t117 = -pkin(4) - t121;
t265 = -t117 / 0.2e1;
t230 = t158 * t160;
t233 = t157 * t163;
t129 = (t230 + t233) * pkin(2);
t264 = t129 / 0.2e1;
t130 = t158 * t254 - t148;
t263 = -t130 / 0.2e1;
t262 = t130 / 0.2e1;
t261 = t148 / 0.2e1;
t150 = -t158 * pkin(3) - pkin(4);
t260 = -t150 / 0.2e1;
t257 = -t159 / 0.2e1;
t256 = t159 / 0.2e1;
t255 = -t162 / 0.2e1;
t63 = t78 * qJD(3);
t81 = t84 * qJD(2);
t252 = -t81 - t63;
t251 = pkin(1) * qJD(1);
t250 = pkin(1) * qJD(2);
t249 = pkin(2) * qJD(2);
t248 = pkin(2) * qJD(3);
t247 = pkin(3) * qJD(3);
t29 = t215 * t79;
t70 = t157 * t120 + t232;
t66 = pkin(8) + t70;
t5 = t66 * t29 + t65 * t78;
t246 = t5 * qJD(1);
t31 = t215 * t85;
t6 = t66 * t31 + t65 * t84;
t245 = t6 * qJD(1);
t244 = t84 * t150;
t243 = t84 * t162;
t149 = t157 * pkin(3) + pkin(8);
t242 = t85 * t149;
t9 = -t69 * t78 + t70 * t79;
t241 = t9 * qJD(1);
t10 = -t69 * t84 + t70 * t85;
t240 = t10 * qJD(1);
t239 = t117 * t162;
t238 = t150 * t162;
t196 = t264 + t271;
t187 = t268 + t196;
t18 = t187 * t162;
t225 = t18 * qJD(1);
t224 = t29 * qJD(1);
t223 = t31 * qJD(1);
t36 = t125 * t131 + t126 * t132;
t222 = t36 * qJD(1);
t221 = t78 * qJD(1);
t220 = t79 * qJD(1);
t219 = t84 * qJD(1);
t218 = t85 * qJD(1);
t217 = t90 * qJD(1);
t91 = t139 / 0.2e1 + (-t253 / 0.2e1 + pkin(2) / 0.2e1) * t163;
t216 = t91 * qJD(1);
t214 = qJD(1) * t159;
t213 = qJD(1) * t162;
t212 = qJD(2) * t159;
t211 = qJD(3) * t159;
t210 = t125 * qJD(1);
t209 = t126 * qJD(1);
t208 = t129 * qJD(2);
t123 = t129 * qJD(3);
t207 = t131 * qJD(1);
t206 = t132 * qJD(1);
t205 = t159 * qJD(5);
t154 = t162 * qJD(5);
t204 = -qJD(1) - qJD(2);
t202 = -t254 / 0.2e1;
t201 = t65 * t214;
t200 = t65 * t213;
t199 = t78 * t214;
t198 = t84 * t214;
t194 = -t230 / 0.2e1;
t191 = pkin(1) * t204;
t190 = pkin(2) * t203;
t83 = t215 * t130;
t189 = t215 * t149;
t122 = pkin(2) * t230 + t157 * t151;
t188 = t202 + t125 / 0.2e1;
t185 = t270 + t260 + t272;
t184 = t267 + t265 + t272;
t182 = t193 * t79;
t181 = t193 * t130;
t180 = t263 + t260 + t265;
t118 = pkin(8) + t122;
t15 = t117 * t129 + t118 * t83;
t174 = t117 * t271 + t65 * t264;
t2 = -t244 / 0.2e1 + t174 + t215 * (t118 * t269 + t66 * t262 - t242 / 0.2e1);
t179 = -t2 * qJD(1) - t15 * qJD(2);
t165 = t121 * t271 + t122 * t270 + t70 * t263 + t69 * t264;
t169 = (t85 * t157 / 0.2e1 + t158 * t268) * pkin(3);
t3 = t169 + t165;
t33 = -t121 * t129 + t122 * t130;
t178 = t3 * qJD(1) - t33 * qJD(2);
t8 = t215 * (t262 + t269 + t267);
t177 = -t8 * qJD(1) - t83 * qJD(2);
t173 = t132 / 0.2e1 + t188;
t24 = pkin(2) * t194 + t266 + t119 / 0.2e1 + t173 * t157;
t176 = -t24 * qJD(1) + t208;
t26 = t261 + (t126 / 0.2e1 - t131 / 0.2e1) * t157 + t173 * t158;
t175 = -t26 * qJD(1) + t130 * qJD(2);
t11 = t184 * t159;
t172 = t11 * qJD(1) - t117 * t212;
t12 = t184 * t162;
t171 = t12 * qJD(1) - qJD(2) * t239;
t17 = t187 * t159;
t170 = -t17 * qJD(1) - t159 * t208;
t20 = t185 * t159;
t45 = t180 * t159;
t168 = t20 * qJD(1) + t45 * qJD(2) - t150 * t211;
t21 = t185 * t162;
t46 = t180 * t162;
t167 = t21 * qJD(1) + t46 * qJD(2) - qJD(3) * t238;
t147 = t159 * t154;
t138 = t146 * qJD(5);
t134 = t238 / 0.2e1;
t133 = t150 * t256;
t128 = t132 * qJD(2);
t127 = t131 * qJD(2);
t124 = t130 * qJD(3);
t116 = t126 * qJD(3);
t115 = t125 * qJD(3);
t110 = t159 * t123;
t94 = t192 * t162 * t159;
t93 = t239 / 0.2e1;
t92 = t117 * t256;
t87 = t202 - t139 / 0.2e1 + (t229 - t226 / 0.2e1) * pkin(1);
t86 = -t90 / 0.2e1 - t183 / 0.2e1 + (-t227 - t228 / 0.2e1) * pkin(1);
t82 = t85 * qJD(2);
t80 = t83 * qJD(3);
t77 = t84 * t212;
t64 = t79 * qJD(3);
t59 = t78 * t211;
t50 = t65 * t162 / 0.2e1;
t49 = t65 * t256;
t48 = t130 * t255 + t134 + t93;
t47 = t130 * t257 + t133 + t92;
t30 = t31 * qJD(2);
t28 = t29 * qJD(3);
t27 = t261 + t236 / 0.2e1 - t231 / 0.2e1 + t235 / 0.2e1 + t188 * t158;
t25 = t237 / 0.2e1 + t266 - t234 / 0.2e1 - t119 / 0.2e1 + (-t233 / 0.2e1 + t194) * pkin(2);
t23 = t79 * t255 + t134 + t50;
t22 = t79 * t257 + t133 + t49;
t19 = -t243 / 0.2e1 - t196 * t162;
t16 = t196 * t159 + t84 * t256;
t14 = t85 * t255 + t50 + t93;
t13 = t85 * t257 + t49 + t92;
t7 = t193 * t85 + t181 + t182;
t4 = t169 - t165;
t1 = t244 / 0.2e1 + t66 * t181 + t118 * t182 + t174 + t215 * t242 / 0.2e1;
t32 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161 * t250, -t164 * t250, 0, 0, 0, 0, 0, 0, 0, 0, -t127 - t116, -t128 + t115, 0, t36 * qJD(2), 0, 0, 0, 0, 0, 0, t252, -t82 - t64, 0, t10 * qJD(2) + t9 * qJD(3), t147, t138, 0, -t147, 0, 0, t252 * t162 + t65 * t205, t65 * t154 + t59 + t77, t30 + t28, t6 * qJD(2) + t5 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161 * t191, t164 * t191, 0, 0, 0, 0, 0, 0, 0, 0, t86 * qJD(3) - t127 - t207, t87 * qJD(3) - t128 - t206, 0, t222 + (-t131 * t163 + t132 * t160) * t249, 0, 0, 0, 0, 0, 0, t25 * qJD(3) - t219 - t81, t27 * qJD(3) - t218 - t82, 0, t240 + (-t84 * t121 + t85 * t122) * qJD(2) + t4 * qJD(3), t147, t138, 0, -t147, 0, 0, t19 * qJD(3) + t13 * qJD(5) + t204 * t243, t16 * qJD(3) + t14 * qJD(5) + t198 + t77, t7 * qJD(3) + t223 + t30, t245 + (t84 * t117 + t118 * t31) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * qJD(2) - t116 - t209, t87 * qJD(2) + t115 + t210, 0, 0, 0, 0, 0, 0, 0, 0, t25 * qJD(2) - t221 - t63, t27 * qJD(2) - t220 - t64, 0, t241 + t4 * qJD(2) + (t157 * t79 - t158 * t78) * t247, t147, t138, 0, -t147, 0, 0, t19 * qJD(2) + t22 * qJD(5) + (-qJD(1) - qJD(3)) * t78 * t162, t16 * qJD(2) + t23 * qJD(5) + t199 + t59, t7 * qJD(2) + t224 + t28, t246 + t1 * qJD(2) + (t78 * t150 + t189 * t79) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t273, t154, -t94, -t205, 0, t13 * qJD(2) + t22 * qJD(3) - t66 * t154 + t201, t14 * qJD(2) + t23 * qJD(3) + t205 * t66 + t200, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161 * t251, t164 * t251, 0, 0, 0, 0, 0, 0, 0, 0, -t90 * qJD(3) + t207, -t91 * qJD(3) + t206, 0, -t222, 0, 0, 0, 0, 0, 0, t24 * qJD(3) + t219, t26 * qJD(3) + t218, 0, -t3 * qJD(3) - t240, t147, t138, 0, -t147, 0, 0, -t18 * qJD(3) - t11 * qJD(5) + t84 * t213, t17 * qJD(3) - t12 * qJD(5) - t198, t8 * qJD(3) - t223, t2 * qJD(3) - t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t248, -t163 * t248, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t124, 0, t33 * qJD(3), t147, t138, 0, -t147, 0, 0, t117 * t205 - t162 * t123, t117 * t154 + t110, t80, t15 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160 * t190 - t217, t163 * t190 - t216, 0, 0, 0, 0, 0, 0, 0, 0, -t123 - t176, -t124 - t175, 0, (-t129 * t158 + t130 * t157) * t247 - t178, t147, t138, 0, -t147, 0, 0, t203 * t162 * t129 + t47 * qJD(5) - t225, t48 * qJD(5) + t110 - t170, -t177 + t80, (t129 * t150 + t130 * t189) * qJD(3) - t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t273, t154, -t94, -t205, 0, t47 * qJD(3) - t118 * t154 - t172, t48 * qJD(3) + t118 * t205 - t171, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 * qJD(2) + t209, t91 * qJD(2) - t210, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * qJD(2) + t221, -t26 * qJD(2) + t220, 0, t3 * qJD(2) - t241, t147, t138, 0, -t147, 0, 0, t18 * qJD(2) - t20 * qJD(5) + t78 * t213, -t17 * qJD(2) - t21 * qJD(5) - t199, -t8 * qJD(2) - t224, -t2 * qJD(2) - t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160 * t249 + t217, t163 * t249 + t216, 0, 0, 0, 0, 0, 0, 0, 0, t176, t175, 0, t178, t147, t138, 0, -t147, 0, 0, -t45 * qJD(5) + t162 * t208 + t225, -t46 * qJD(5) + t170, t177, t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t138, 0, -t147, 0, 0, t150 * t205, t150 * t154, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t273, t154, -t94, -t205, 0, -t149 * t154 - t168, t149 * t205 - t167, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, -t154, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t273, 0, t94, 0, 0, t11 * qJD(2) + t20 * qJD(3) - t201, t12 * qJD(2) + t21 * qJD(3) - t200, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t273, 0, t94, 0, 0, t45 * qJD(3) + t172, t46 * qJD(3) + t171, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t273, 0, t94, 0, 0, t168, t167, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t32;
