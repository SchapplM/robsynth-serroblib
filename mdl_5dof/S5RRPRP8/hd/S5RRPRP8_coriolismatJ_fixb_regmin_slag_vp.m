% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:24
% EndTime: 2021-01-15 20:52:32
% DurationCPUTime: 2.06s
% Computational Cost: add. (1989->195), mult. (3530->241), div. (0->0), fcn. (3464->4), ass. (0->167)
t149 = sin(qJ(2));
t150 = cos(qJ(2));
t239 = sin(qJ(4));
t240 = cos(qJ(4));
t114 = -t149 * t239 - t150 * t240;
t246 = pkin(6) - pkin(7);
t128 = t246 * t149;
t129 = t246 * t150;
t51 = t239 * t128 + t240 * t129;
t19 = t114 * qJ(5) + t51;
t217 = t19 * qJD(4);
t263 = t19 * qJD(2) - t217;
t187 = -t239 / 0.2e1;
t115 = t149 * t240 - t150 * t239;
t100 = t115 * qJ(5);
t117 = t239 * t129;
t175 = t240 * t128;
t207 = t117 - t175;
t44 = t100 + t207;
t261 = t44 * t187;
t189 = qJD(2) - qJD(4);
t107 = t114 ^ 2;
t250 = t115 ^ 2;
t48 = t250 - t107;
t260 = t48 * qJD(1);
t247 = -pkin(2) - pkin(3);
t122 = t239 * qJ(3) - t240 * t247;
t245 = t175 / 0.2e1;
t66 = t245 - t175 / 0.2e1;
t259 = t66 * qJD(1) - t122 * qJD(2);
t258 = t189 * t51;
t167 = t240 * qJD(3);
t257 = -t167 - t259;
t108 = t122 * qJD(4);
t256 = t108 + t259;
t190 = t150 * qJD(3);
t206 = -t150 * pkin(2) - t149 * qJ(3);
t255 = t206 * qJD(2) + t190;
t59 = t189 * t114;
t58 = t189 * t115;
t252 = t66 * qJD(2);
t251 = t66 * qJD(4);
t249 = -pkin(4) / 0.2e1;
t248 = t19 / 0.2e1;
t121 = -pkin(4) - t122;
t244 = -t121 / 0.2e1;
t243 = t121 / 0.2e1;
t241 = t149 / 0.2e1;
t238 = t114 * pkin(4);
t237 = t115 * pkin(4);
t180 = t19 * t240;
t181 = t44 * t239;
t236 = t181 / 0.2e1 + t180 / 0.2e1;
t124 = -pkin(1) + t206;
t99 = t150 * pkin(3) - t124;
t60 = t99 - t238;
t144 = t150 * qJ(3);
t105 = t247 * t149 + t144;
t61 = t105 - t237;
t3 = t60 * t61;
t234 = t3 * qJD(1);
t231 = t60 * t114;
t230 = t60 * t115;
t170 = t122 / 0.2e1 + t243;
t12 = (t249 - t170) * t114;
t228 = t12 * qJD(1);
t123 = t240 * qJ(3) + t239 * t247;
t227 = t123 * t114;
t13 = t114 * t19 + t115 * t44;
t226 = t13 * qJD(1);
t14 = t61 * t115 - t231;
t225 = t14 * qJD(1);
t15 = -t61 * t114 - t230;
t224 = t15 * qJD(1);
t154 = (-pkin(2) / 0.2e1 - pkin(3) / 0.2e1) * t149 + t144 / 0.2e1;
t16 = -t227 / 0.2e1 + (t249 + t243) * t115 + t154;
t223 = t16 * qJD(1);
t20 = -t250 * pkin(4) - t231;
t222 = t20 * qJD(1);
t21 = t114 * t237 - t230;
t221 = t21 * qJD(1);
t26 = t105 * t115 - t99 * t114;
t219 = t26 * qJD(1);
t27 = -t105 * t114 - t99 * t115;
t218 = t27 * qJD(1);
t151 = t240 * t115 / 0.2e1 + t114 * t187;
t52 = t241 + t151;
t214 = t52 * qJD(1);
t62 = t107 + t250;
t213 = t62 * qJD(1);
t127 = t149 * pkin(2) - t144;
t63 = t124 * t150 + t127 * t149;
t212 = t63 * qJD(1);
t64 = -t124 * t149 + t127 * t150;
t211 = t64 * qJD(1);
t112 = t123 * qJD(4);
t194 = t123 * qJD(2);
t208 = t194 - t112;
t205 = qJD(1) * t149;
t204 = qJD(1) * t150;
t203 = qJD(2) * qJ(3);
t202 = qJD(3) * t149;
t201 = t114 * qJD(1);
t200 = t114 * qJD(4);
t199 = t114 * qJD(5);
t198 = t115 * qJD(1);
t197 = t115 * qJD(4);
t196 = t115 * qJD(5);
t148 = t149 ^ 2;
t130 = t150 ^ 2 - t148;
t193 = t130 * qJD(1);
t192 = t148 * qJD(1);
t191 = t149 * qJD(2);
t140 = t150 * qJD(2);
t188 = pkin(4) * t248;
t186 = pkin(1) * t205;
t185 = pkin(1) * t204;
t184 = pkin(4) * t198;
t183 = pkin(6) * t191;
t182 = pkin(6) * t140;
t179 = t60 * t205;
t178 = t99 * t201;
t177 = t99 * t198;
t176 = t121 * t239;
t173 = t114 * t198;
t172 = t124 * t127 * qJD(1);
t171 = t124 * t205;
t169 = qJD(4) * t239;
t168 = t240 * qJD(2);
t166 = t239 * qJD(2);
t165 = t239 * qJD(3);
t50 = -t117 + 0.2e1 * t245;
t163 = pkin(4) * t187;
t162 = -t180 / 0.2e1;
t159 = -t194 - t165;
t6 = pkin(4) * t230;
t8 = t162 + t236 + t261;
t158 = t6 * qJD(1) + t8 * qJD(3);
t1 = t170 * t19 + t188;
t33 = (-t121 - t122) * t123;
t157 = t1 * qJD(1) - t33 * qJD(2);
t152 = -t176 / 0.2e1 + t122 * t187;
t22 = t163 - t152;
t156 = -t8 * qJD(1) + t22 * qJD(2);
t57 = t123 * t240 - t176;
t153 = t240 * t248 - t261;
t9 = t162 - t181 / 0.2e1 + t153;
t155 = t9 * qJD(1) - t57 * qJD(2);
t131 = t149 * t204;
t126 = -qJD(4) * t240 + t168;
t125 = t166 - t169;
t84 = t149 * t201;
t83 = t114 * t202;
t82 = t149 * t198;
t81 = t115 * t202;
t74 = t167 - t108;
t73 = t165 + t112;
t53 = t241 - t151;
t32 = t240 * t114 + t239 * t115;
t23 = t163 + t152;
t18 = -t100 + t50;
t17 = t227 / 0.2e1 + t115 * t244 - t237 / 0.2e1 + t154;
t11 = t238 / 0.2e1 - t170 * t114;
t10 = t153 + t236;
t7 = t8 * qJD(4);
t2 = t188 + (-t122 / 0.2e1 + t244) * t19;
t4 = [0, 0, 0, t149 * t140, t130 * qJD(2), 0, 0, 0, -pkin(1) * t191, -pkin(1) * t140, -t64 * qJD(2) + t149 * t190, 0, -t63 * qJD(2) + t148 * qJD(3), (qJD(2) * t127 - t202) * t124, -t114 * t58, t189 * t48, 0, 0, 0, t27 * qJD(2) + t99 * t197 - t83, t26 * qJD(2) + t99 * t200 + t81, t15 * qJD(2) - t21 * qJD(4) - t83, t14 * qJD(2) - t20 * qJD(4) + t81, t62 * qJD(5), t3 * qJD(2) + t6 * qJD(4) + t13 * qJD(5) + t60 * t202; 0, 0, 0, t131, t193, t140, -t191, 0, -t182 - t186, t183 - t185, -t182 - t211, t255, -t183 - t212, t255 * pkin(6) + t172, -t173, t260, t59, -t58, 0, t218 - t258, t207 * qJD(2) + t50 * qJD(4) + t219, t224 - t263, qJD(2) * t44 + t18 * qJD(4) + t225, (t121 * t114 + t123 * t115) * qJD(2) + t32 * qJD(3) + t11 * qJD(4), t234 + (t121 * t19 + t44 * t123) * qJD(2) + t10 * qJD(3) + t2 * qJD(4) + t17 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t140, t192, -t171 + t182, 0, 0, 0, 0, 0, -t84, t82, -t84, t82, t32 * qJD(2), t10 * qJD(2) + t53 * qJD(5) + t179 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, -t260, -t59, t58, 0, t177 + t258, t50 * qJD(2) + t207 * qJD(4) + t178, -t221 + t263, t18 * qJD(2) + t44 * qJD(4) - t222, -pkin(4) * t200 + t11 * qJD(2), -pkin(4) * t217 + t2 * qJD(2) + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t17 * qJD(2) + t53 * qJD(3) + t226; 0, 0, 0, -t131, -t193, 0, 0, 0, t186, t185, t211, 0, t212, -t172, t173, -t260, 0, 0, 0, -t218, -t219 + t251, t196 - t224, t199 - t225 + t251, t12 * qJD(4), -t9 * qJD(3) - t1 * qJD(4) - t16 * qJD(5) - t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, 0, 0, 0, 0, t73, t74, t73, t74, 0, t57 * qJD(3) + t33 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t203, 0, 0, 0, 0, 0, t166, t168, t166, t168, 0, t23 * qJD(4) - t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t256, t208, t256, t228, -pkin(4) * t112 + t23 * qJD(3) - t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t201, 0, -t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, 0, -t192, t171, 0, 0, 0, 0, 0, t84, -t82, t84, -t82, 0, t9 * qJD(2) - t52 * qJD(5) - t179 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t203, 0, 0, 0, 0, 0, -t125, -t126, -t125, -t126, 0, -t22 * qJD(4) + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t126, t125, t126, 0, -pkin(4) * t169 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t260, 0, 0, 0, -t177, -t178 - t252, -t196 + t221, -t199 + t222 - t252, -t12 * qJD(2), -pkin(4) * t196 + t1 * qJD(2) - t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t257, t159, t257, -t228, t22 * qJD(3) + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t168, -t166, -t168, 0, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t201, 0, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t59, -t213, pkin(4) * t197 + t16 * qJD(2) + t52 * qJD(3) - t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t201, 0, t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t201, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
