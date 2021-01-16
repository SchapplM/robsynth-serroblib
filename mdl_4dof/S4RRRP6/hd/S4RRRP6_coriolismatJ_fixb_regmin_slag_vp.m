% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:17
% EndTime: 2021-01-15 14:39:24
% DurationCPUTime: 1.80s
% Computational Cost: add. (1340->227), mult. (3026->344), div. (0->0), fcn. (2467->4), ass. (0->196)
t129 = cos(qJ(2));
t127 = sin(qJ(2));
t128 = cos(qJ(3));
t215 = t128 * qJ(4);
t147 = -t129 * pkin(2) - t127 * pkin(6);
t142 = -pkin(1) + t147;
t77 = t128 * t142;
t148 = -t127 * t215 + t77;
t126 = sin(qJ(3));
t243 = pkin(5) * t126;
t38 = (-pkin(3) - t243) * t129 + t148;
t217 = t126 * t129;
t184 = pkin(5) * t217;
t44 = -t148 + t184;
t251 = t38 + t44;
t122 = t126 ^ 2;
t124 = t128 ^ 2;
t154 = t122 / 0.2e1 - t124 / 0.2e1;
t110 = t124 - t122;
t193 = t127 * qJD(1);
t162 = t128 * t193;
t250 = t110 * qJD(2) - 0.2e1 * t126 * t162;
t218 = t126 * t127;
t116 = pkin(5) * t218;
t212 = t129 * qJ(4);
t240 = t127 * pkin(3);
t237 = t129 * pkin(6);
t241 = t127 * pkin(2);
t95 = -t237 + t241;
t87 = t128 * t95;
t41 = -t128 * t212 + t116 + t240 + t87;
t249 = t41 / 0.2e1;
t242 = t126 * pkin(3);
t171 = pkin(5) + t242;
t88 = t171 * t127;
t248 = -t88 / 0.2e1;
t236 = qJ(4) + pkin(6);
t92 = t236 * t128;
t247 = -t92 / 0.2e1;
t244 = t127 / 0.2e1;
t239 = t128 * pkin(3);
t238 = t129 * pkin(5);
t174 = t44 / 0.2e1 + t38 / 0.2e1;
t175 = -t129 * pkin(3) / 0.2e1;
t3 = (t175 + t174) * t128;
t235 = t3 * qJD(1);
t234 = t38 * t128;
t213 = t128 * t129;
t182 = pkin(5) * t213;
t51 = t126 * t142 + t182;
t45 = -qJ(4) * t218 + t51;
t233 = t45 * t129;
t216 = t127 * t128;
t183 = pkin(5) * t216;
t86 = t126 * t95;
t46 = -t126 * t212 - t183 + t86;
t89 = t171 * t129;
t5 = t38 * t41 + t45 * t46 + t88 * t89;
t232 = t5 * qJD(1);
t6 = (t127 * t41 + t129 * t38) * t128 + (t127 * t46 + t233) * t126;
t231 = t6 * qJD(1);
t7 = t251 * t218;
t230 = t7 * qJD(1);
t177 = t88 * t216;
t8 = pkin(3) * t177 - t251 * t45;
t229 = t8 * qJD(1);
t228 = t88 * t128;
t227 = t89 * t128;
t91 = t236 * t126;
t226 = t91 * t127;
t225 = t91 * t129;
t224 = t92 * t127;
t223 = t92 * t129;
t11 = -t38 * t127 + t41 * t129 - t88 * t217 - t89 * t218;
t222 = t11 * qJD(1);
t119 = -pkin(2) - t239;
t221 = t119 * t126;
t220 = t119 * t128;
t12 = (t46 + t228) * t129 + (-t45 + t227) * t127;
t219 = t12 * qJD(1);
t123 = t127 ^ 2;
t214 = t128 * t123;
t17 = (t45 * t126 + t234) * t127;
t211 = t17 * qJD(1);
t173 = t126 * t214;
t18 = -pkin(3) * t173 - t177 - t233;
t210 = t18 * qJD(1);
t50 = -t77 + t184;
t19 = t50 * t127 + (-t116 + t87) * t129;
t209 = t19 * qJD(1);
t20 = t86 * t129 + (-t51 + t182) * t127;
t208 = t20 * qJD(1);
t21 = t124 * t123 * pkin(3) - t44 * t129 - t88 * t218;
t207 = t21 * qJD(1);
t39 = -t123 * t243 - t50 * t129;
t206 = t39 * qJD(1);
t40 = -pkin(5) * t214 - t51 * t129;
t205 = t40 * qJD(1);
t204 = t45 * qJD(3);
t109 = t124 + t122;
t82 = t109 * t123;
t203 = t82 * qJD(1);
t125 = t129 ^ 2;
t111 = t125 - t123;
t84 = t111 * t126;
t202 = t84 * qJD(1);
t85 = t128 * t125 - t214;
t201 = t85 * qJD(1);
t200 = t92 * qJD(3);
t199 = qJD(3) * t126;
t121 = qJD(3) * t128;
t198 = qJD(3) * t129;
t197 = t109 * qJD(2);
t196 = t111 * qJD(1);
t195 = t126 * qJD(2);
t194 = t126 * qJD(4);
t192 = t127 * qJD(2);
t191 = t128 * qJD(2);
t190 = t128 * qJD(4);
t189 = t129 * qJD(1);
t188 = t129 * qJD(2);
t187 = t129 * qJD(4);
t75 = t87 / 0.2e1;
t186 = t75 + t116 / 0.2e1;
t185 = t126 * t239;
t181 = pkin(1) * t193;
t180 = pkin(1) * t189;
t179 = pkin(3) * t199;
t178 = pkin(3) * t121;
t176 = t242 / 0.2e1;
t172 = t126 * t248;
t170 = t126 * t198;
t169 = t128 * t198;
t168 = t126 * t121;
t167 = t126 * t191;
t166 = t126 * t187;
t165 = t127 * t188;
t164 = t127 * t189;
t163 = t127 * t191;
t161 = t127 * t190;
t160 = t128 * t189;
t159 = t128 * t187;
t158 = -t220 / 0.2e1;
t157 = t216 / 0.2e1;
t156 = -t213 / 0.2e1;
t155 = t212 / 0.2e1;
t153 = -t127 * t185 + t225 / 0.2e1;
t152 = -qJD(3) + t189;
t150 = t126 * t163;
t149 = -t224 / 0.2e1 - t38 / 0.2e1;
t135 = t127 * t158 + t172;
t1 = t174 * t92 + (t249 + t135) * pkin(3);
t22 = pkin(3) * t221;
t146 = -t1 * qJD(1) + t22 * qJD(2);
t139 = (t45 / 0.2e1 + t226 / 0.2e1) * t128;
t10 = -t238 / 0.2e1 + t139 + (t175 + t149) * t126;
t48 = t91 * t126 + t92 * t128;
t145 = t10 * qJD(1) + t48 * qJD(2);
t13 = t172 + (-t215 / 0.2e1 + t247) * t129 + (t158 + (0.1e1 - t154) * pkin(3)) * t127 + t186;
t64 = t185 - t221;
t144 = -t13 * qJD(1) - t64 * qJD(2);
t74 = -t86 / 0.2e1;
t15 = t74 + (pkin(5) * t244 + t248) * t128 + (t119 * t244 + t155) * t126 + t153;
t72 = t122 * pkin(3) + t220;
t143 = -t15 * qJD(1) + t72 * qJD(2);
t141 = t152 * t127;
t140 = t237 / 0.2e1 - t241 / 0.2e1;
t133 = t140 * t128;
t43 = -t87 / 0.2e1 + t133;
t138 = pkin(2) * t195 - t43 * qJD(1);
t134 = t140 * t126;
t42 = t86 / 0.2e1 - t134;
t137 = pkin(2) * t191 - t42 * qJD(1);
t71 = t154 * t127;
t136 = -t71 * qJD(1) + t167;
t55 = t128 * t141;
t81 = t162 + t195;
t80 = t126 * t193 - t191;
t132 = qJD(1) * t173 + t71 * qJD(2);
t83 = t110 * t123;
t131 = t83 * qJD(1) + 0.2e1 * t150;
t117 = t192 / 0.2e1;
t76 = (t189 - qJD(3) / 0.2e1) * t127;
t68 = t81 * pkin(3);
t67 = t71 * qJD(3);
t54 = t81 * t129;
t53 = t80 * t129;
t52 = t126 * t141;
t31 = t116 + t75 + t133;
t30 = -t134 + t74 + t183;
t16 = -t119 * t218 / 0.2e1 + t228 / 0.2e1 + pkin(5) * t157 + t74 + t126 * t155 - t153;
t14 = t223 / 0.2e1 + qJ(4) * t156 - t135 + t186 + (0.1e1 + t154) * t240;
t9 = t238 / 0.2e1 + t129 * t176 + t139 + t149 * t126;
t4 = -t44 * t128 / 0.2e1 - t234 / 0.2e1 + pkin(3) * t156;
t2 = t88 * t176 + t251 * t247 + (t119 * t157 + t249) * pkin(3);
t23 = [0, 0, 0, t165, t111 * qJD(2), 0, 0, 0, -pkin(1) * t192, -pkin(1) * t188, -t123 * t168 + t124 * t165, -t83 * qJD(3) - 0.2e1 * t129 * t150, -t85 * qJD(2) + t127 * t170, t84 * qJD(2) + t127 * t169, -t165, -t19 * qJD(2) - t40 * qJD(3), t20 * qJD(2) + t39 * qJD(3), -t11 * qJD(2) - t18 * qJD(3) + t127 * t159, t12 * qJD(2) + t21 * qJD(3) - t127 * t166, -t6 * qJD(2) + t7 * qJD(3) + t82 * qJD(4), t5 * qJD(2) + t8 * qJD(3) - t17 * qJD(4); 0, 0, 0, t164, t196, t188, -t192, 0, -pkin(5) * t188 - t181, pkin(5) * t192 - t180, -t67 + (t124 * t193 + t167) * t129, -0.2e1 * t127 * t168 + t250 * t129, t126 * t192 - t201, t163 + t202, -t76, -t209 + (t126 * t147 - t182) * qJD(2) + t31 * qJD(3), t208 + (t128 * t147 + t184) * qJD(2) + t30 * qJD(3), -t222 + (t119 * t217 - t226 - t227) * qJD(2) + t14 * qJD(3) + t166, t219 + (t119 * t213 + t89 * t126 - t224) * qJD(2) + t16 * qJD(3) + t159, -t231 + ((t46 + t225) * t128 + (-t41 - t223) * t126) * qJD(2) + t4 * qJD(3), t232 + (t89 * t119 - t41 * t91 + t46 * t92) * qJD(2) + t2 * qJD(3) + t9 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131, t52, t55, t117, t31 * qJD(2) - t51 * qJD(3) - t205, t30 * qJD(2) + t50 * qJD(3) + t206, t14 * qJD(2) - t204 - t210, t16 * qJD(2) + t44 * qJD(3) + t207, t4 * qJD(2) + t127 * t179 + t230, -pkin(3) * t204 + t2 * qJD(2) + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t53, t203, t9 * qJD(2) - t211; 0, 0, 0, -t164, -t196, 0, 0, 0, t181, t180, -t124 * t164 - t67, 0.2e1 * t126 * t55, -t169 + t201, t170 - t202, t76, t43 * qJD(3) + t209, t42 * qJD(3) - t208, -t13 * qJD(3) + t222, -t15 * qJD(3) - t219, -qJD(3) * t3 + t231, -qJD(3) * t1 + qJD(4) * t10 - t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t110 * qJD(3), 0, 0, 0, -pkin(2) * t199, -pkin(2) * t121, -t64 * qJD(3), t72 * qJD(3), t109 * qJD(4), t22 * qJD(3) + t48 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t250, t121 - t160, t152 * t126, -t193 / 0.2e1, -pkin(6) * t121 - t138, pkin(6) * t199 - t137, t144 - t200, t91 * qJD(3) + t143, -t178 - t235, -pkin(3) * t200 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t131, -t53, -t54, t117, -t43 * qJD(2) + t205, -t42 * qJD(2) - t206, t13 * qJD(2) - t161 + t210, t15 * qJD(2) + t127 * t194 - t207, qJD(2) * t3 - t230, -pkin(3) * t161 + t1 * qJD(2) - t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t250, t160, -t126 * t189, t193 / 0.2e1, t138, t137, -t144 - t194, -t143 - t190, t235, -pkin(3) * t194 - t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t80, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t52, -t203, -t10 * qJD(2) + t127 * t178 + t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t121, -t197, -t145 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t80, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t23;
