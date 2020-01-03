% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:16
% EndTime: 2019-12-31 17:26:21
% DurationCPUTime: 1.73s
% Computational Cost: add. (1792->191), mult. (4010->291), div. (0->0), fcn. (4022->6), ass. (0->177)
t203 = qJD(2) + qJD(3);
t149 = sin(qJ(2));
t151 = cos(qJ(2));
t249 = sin(qJ(3));
t250 = cos(qJ(3));
t121 = t249 * t149 - t250 * t151;
t139 = t249 * t151;
t140 = t250 * t149;
t123 = -t140 - t139;
t262 = t203 * t123;
t173 = t121 * t262;
t150 = cos(qJ(4));
t144 = -t151 * pkin(2) - pkin(1);
t245 = t123 * pkin(7);
t248 = t121 * pkin(3);
t176 = t245 + t248;
t159 = t144 + t176;
t148 = sin(qJ(4));
t256 = pkin(5) + pkin(6);
t129 = t256 * t149;
t130 = t256 * t151;
t160 = -t249 * t129 + t250 * t130;
t264 = t160 * t148;
t25 = -t150 * t159 + t264;
t269 = (t25 - t264) * t123;
t263 = t160 * t150;
t26 = t148 * t159 + t263;
t268 = (t26 - t263) * t123;
t175 = t250 * t129 + t249 * t130;
t267 = t203 * t175;
t266 = t203 * t160;
t265 = 0.2e1 * t123;
t73 = t203 * t121;
t120 = t123 ^ 2;
t202 = -t121 ^ 2 + t120;
t199 = t249 * pkin(2);
t142 = t199 + pkin(7);
t227 = t123 * t142;
t200 = t250 * pkin(2);
t143 = -t200 - pkin(3);
t228 = t121 * t143;
t260 = t245 / 0.2e1 + t248 / 0.2e1 + t227 / 0.2e1 - t228 / 0.2e1 + (-t250 * t121 / 0.2e1 - t249 * t123 / 0.2e1) * pkin(2);
t226 = t148 * t150;
t192 = qJD(1) * t226;
t147 = t150 ^ 2;
t257 = t148 ^ 2;
t52 = (-t257 / 0.2e1 + t147 / 0.2e1) * t123;
t259 = t120 * t192 + t203 * t52;
t201 = -t147 + t257;
t40 = -0.2e1 * t123 * t192 + t201 * t203;
t255 = t142 / 0.2e1;
t254 = t143 / 0.2e1;
t253 = -t148 / 0.2e1;
t252 = -t150 / 0.2e1;
t251 = t150 / 0.2e1;
t247 = t121 * pkin(7);
t246 = t123 * pkin(3);
t244 = t149 * pkin(2);
t242 = pkin(3) * qJD(3);
t75 = -t246 + t247;
t59 = t75 + t244;
t238 = t150 * t59;
t1 = t238 * t121 + t269;
t241 = t1 * qJD(1);
t240 = t148 * t59;
t2 = -t240 * t121 + t268;
t236 = t2 * qJD(1);
t56 = t150 * t121;
t3 = t75 * t56 + t269;
t235 = t3 * qJD(1);
t53 = t148 * t121;
t4 = -t75 * t53 + t268;
t234 = t4 * qJD(1);
t179 = -t200 / 0.2e1;
t169 = t179 - t143 / 0.2e1;
t152 = (t255 - t199 / 0.2e1 - pkin(7) / 0.2e1) * t123 + (-pkin(3) / 0.2e1 + t169) * t121;
t5 = t152 * t148;
t233 = t5 * qJD(1);
t230 = t175 * t148;
t229 = t175 * t150;
t18 = t25 * t121 + t123 * t230;
t225 = t18 * qJD(1);
t19 = -t26 * t121 - t123 * t229;
t224 = t19 * qJD(1);
t31 = t202 * t148;
t223 = t31 * qJD(1);
t32 = t202 * t150;
t222 = t32 * qJD(1);
t221 = t202 * qJD(1);
t47 = t121 * t244 - t144 * t123;
t220 = t47 * qJD(1);
t48 = -t144 * t121 - t123 * t244;
t219 = t48 * qJD(1);
t218 = t52 * qJD(1);
t217 = t53 * qJD(1);
t51 = t56 * qJD(1);
t62 = t201 * t120;
t216 = t62 * qJD(1);
t215 = qJD(1) * t121;
t214 = qJD(1) * t123;
t213 = qJD(1) * t144;
t212 = qJD(1) * t151;
t211 = qJD(2) * t143;
t210 = qJD(3) * t144;
t209 = qJD(4) * t148;
t208 = qJD(4) * t150;
t115 = t140 / 0.2e1 + t139 / 0.2e1;
t207 = t115 * qJD(1);
t135 = -t149 ^ 2 + t151 ^ 2;
t206 = t135 * qJD(1);
t205 = t149 * qJD(2);
t204 = t151 * qJD(2);
t198 = pkin(1) * t149 * qJD(1);
t197 = pkin(1) * t212;
t191 = qJD(4) * t121 * t123;
t76 = t121 * t214;
t190 = t121 * t213;
t189 = t123 * t213;
t138 = t148 * t208;
t188 = t149 * t212;
t187 = t250 * qJD(2);
t186 = t250 * qJD(3);
t185 = t249 * qJD(2);
t184 = t249 * qJD(3);
t182 = t203 * t150;
t181 = pkin(2) * t184;
t180 = pkin(2) * t185;
t174 = t148 * t182;
t172 = t203 * t226;
t171 = t227 - t228;
t170 = t123 * (qJD(4) + t215);
t168 = t247 / 0.2e1 - t246 / 0.2e1;
t166 = t121 * t255 + t123 * t254;
t158 = t59 / 0.2e1 + t166;
t9 = t158 * t148;
t167 = -t9 * qJD(1) - t150 * t211;
t11 = t158 * t150;
t165 = t11 * qJD(1) - t148 * t211;
t164 = t150 * t170;
t34 = -t115 * qJD(4) + t76;
t163 = pkin(3) / 0.2e1 + t169;
t162 = t75 / 0.2e1 + t168;
t161 = t174 * t265;
t8 = t152 * t150;
t157 = -t8 * qJD(1) - t148 * t180;
t13 = t162 * t148;
t92 = t163 * t150;
t156 = -t13 * qJD(1) + t92 * qJD(2) + t150 * t242;
t15 = t162 * t150;
t91 = t163 * t148;
t155 = t15 * qJD(1) + t91 * qJD(2) + t148 * t242;
t154 = (-t185 - t184) * pkin(2);
t133 = t148 * t181;
t128 = t201 * qJD(4);
t94 = pkin(3) * t252 + t143 * t251 + t150 * t179;
t93 = pkin(3) * t253 + (t179 + t254) * t148;
t60 = t203 * t115;
t49 = t52 * qJD(4);
t38 = t208 + t51;
t37 = -t209 - t217;
t30 = t172 - t218;
t29 = -t174 + t218;
t28 = 0.2e1 * t148 * t164;
t27 = -t147 * t76 - t49;
t24 = t56 * qJD(4) - t222;
t23 = -t53 * qJD(4) + t223;
t22 = -t49 + (t147 * t214 - t174) * t121;
t21 = -t148 * t262 + t222;
t20 = -t123 * t182 - t223;
t17 = t201 * t73 + (qJD(4) - t215) * t226 * t265;
t16 = -t168 * t150 + t75 * t251 + t230;
t14 = t168 * t148 + t75 * t253 + t229;
t12 = t230 / 0.2e1 - t175 * t253 + t238 / 0.2e1 - t166 * t150;
t10 = t229 / 0.2e1 - t175 * t252 - t240 / 0.2e1 + t166 * t148;
t7 = t260 * t150 + t264;
t6 = t260 * t148 - t263;
t33 = [0, 0, 0, t149 * t204, t135 * qJD(2), 0, 0, 0, -pkin(1) * t205, -pkin(1) * t204, t173, -t203 * t202, 0, 0, 0, t47 * qJD(2) - t123 * t210, t48 * qJD(2) - t121 * t210, -t120 * t138 + t147 * t173, t62 * qJD(4) - t121 * t161, t148 * t191 + t203 * t32, t150 * t191 - t203 * t31, -t173, t1 * qJD(2) + t3 * qJD(3) + t19 * qJD(4), t2 * qJD(2) + t4 * qJD(3) + t18 * qJD(4); 0, 0, 0, t188, t206, t204, -t205, 0, -pkin(5) * t204 - t198, pkin(5) * t205 - t197, t76, -t221, -t73, t262, 0, t220 - t266, t219 + t267, t22, t17, t21, t20, -t34, t241 + (t148 * t171 - t263) * qJD(2) + t6 * qJD(3) + t12 * qJD(4), t236 + (t150 * t171 + t264) * qJD(2) + t7 * qJD(3) + t10 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t221, -t73, t262, 0, -t189 - t266, -t190 + t267, t22, t17, t21, t20, -t34, t235 + t6 * qJD(2) + (t148 * t176 - t263) * qJD(3) + t16 * qJD(4), t234 + t7 * qJD(2) + (t150 * t176 + t264) * qJD(3) + t14 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, t172 * t265 + t216, t148 * t170, t164, t60, t12 * qJD(2) + t16 * qJD(3) - t26 * qJD(4) + t224, t10 * qJD(2) + t14 * qJD(3) + t25 * qJD(4) + t225; 0, 0, 0, -t188, -t206, 0, 0, 0, t198, t197, -t76, t221, 0, 0, 0, -t220, -t219, t27, t28, t24, t23, t34, t5 * qJD(3) - t11 * qJD(4) - t241, t8 * qJD(3) + t9 * qJD(4) - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -pkin(2) * t186, t138, -t128, 0, 0, 0, t143 * t209 - t150 * t181, t143 * t208 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, (-t187 - t186) * pkin(2), t138, -t128, 0, 0, 0, t93 * qJD(4) + t150 * t154 + t233, t94 * qJD(4) + t133 - t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t40, t38, t37, -t207, t93 * qJD(3) - t142 * t208 - t165, t94 * qJD(3) + t142 * t209 - t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t221, 0, 0, 0, t189, t190, t27, t28, t24, t23, t34, -t5 * qJD(2) - t15 * qJD(4) - t235, -t8 * qJD(2) + t13 * qJD(4) - t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, pkin(2) * t187, t138, -t128, 0, 0, 0, -t91 * qJD(4) + t150 * t180 - t233, -t92 * qJD(4) + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, -t128, 0, 0, 0, -pkin(3) * t209, -pkin(3) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t40, t38, t37, -t207, -pkin(7) * t208 - t155, pkin(7) * t209 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, -t161 - t216, -t148 * t76 - t203 * t56, -t150 * t76 + t203 * t53, t60, t11 * qJD(2) + t15 * qJD(3) - t224, -t9 * qJD(2) - t13 * qJD(3) - t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t40, -t51, t217, t207, t91 * qJD(3) + t165, t92 * qJD(3) + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t40, -t51, t217, t207, t155, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t33;
