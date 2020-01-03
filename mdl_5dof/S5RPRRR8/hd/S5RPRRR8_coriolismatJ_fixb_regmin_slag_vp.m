% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:08
% EndTime: 2019-12-31 19:06:13
% DurationCPUTime: 1.50s
% Computational Cost: add. (1534->192), mult. (2935->249), div. (0->0), fcn. (2856->6), ass. (0->152)
t178 = qJD(1) - qJD(3);
t177 = qJD(4) + qJD(5);
t143 = cos(qJ(5));
t144 = cos(qJ(4));
t210 = t143 * t144;
t140 = sin(qJ(5));
t141 = sin(qJ(4));
t214 = t140 * t141;
t115 = -t210 + t214;
t211 = t143 * t141;
t213 = t140 * t144;
t117 = t211 + t213;
t166 = t115 ^ 2 - t117 ^ 2;
t243 = t177 * t166;
t11 = t178 * t166;
t242 = t177 * t115;
t241 = t177 * t117;
t142 = sin(qJ(3));
t145 = cos(qJ(3));
t234 = -pkin(1) - pkin(2);
t121 = t145 * qJ(2) + t142 * t234;
t164 = t178 * t121;
t132 = -t141 ^ 2 + t144 ^ 2;
t239 = t178 * t132;
t152 = t213 / 0.2e1 + t211 / 0.2e1;
t148 = t152 * t145;
t208 = t145 * t117;
t56 = t208 / 0.2e1 - t148;
t119 = -pkin(7) + t121;
t228 = pkin(8) - t119;
t81 = t228 * t141;
t82 = t228 * t144;
t238 = t56 * qJD(2) + t177 * (-t140 * t81 + t143 * t82);
t151 = -t210 / 0.2e1 + t214 / 0.2e1;
t146 = t151 * t145;
t209 = t145 * t115;
t53 = -t209 / 0.2e1 + t146;
t237 = t53 * qJD(2) + t177 * (-t140 * t82 - t143 * t81);
t233 = pkin(7) + pkin(8);
t130 = t233 * t141;
t131 = t233 * t144;
t49 = (t117 / 0.2e1 - t152) * t145;
t197 = t49 * qJD(2);
t236 = -t197 + t177 * (t140 * t130 - t143 * t131);
t52 = (-t115 / 0.2e1 + t151) * t145;
t196 = t52 * qJD(2);
t235 = -t196 + t177 * (t143 * t130 + t140 * t131);
t232 = t144 / 0.2e1;
t231 = pkin(3) * t141;
t230 = pkin(4) * t141;
t229 = pkin(4) * t144;
t136 = -pkin(3) - t229;
t216 = t136 * t115;
t120 = qJ(2) * t142 - t145 * t234;
t118 = pkin(3) + t120;
t105 = t118 + t229;
t218 = t105 * t115;
t227 = -t218 / 0.2e1 + t216 / 0.2e1;
t168 = t120 * t232;
t212 = t141 * t120;
t226 = -t140 * t212 / 0.2e1 + t143 * t168;
t225 = pkin(4) * qJD(5);
t224 = qJD(4) * pkin(4);
t50 = (-t117 / 0.2e1 - t152) * t145;
t223 = qJD(1) * t50;
t51 = (t115 / 0.2e1 + t151) * t145;
t222 = qJD(1) * t51;
t57 = t142 * t115;
t221 = qJD(1) * t57;
t59 = t142 * t117;
t220 = qJD(1) * t59;
t26 = t117 * t115;
t219 = qJD(3) * t26;
t217 = t105 * t117;
t215 = t136 * t117;
t176 = t115 * t230;
t29 = -t176 + t217;
t203 = t29 * qJD(1);
t175 = t117 * t230;
t30 = -t175 - t218;
t202 = t30 * qJD(1);
t195 = qJ(2) * qJD(1);
t194 = qJD(1) * t117;
t193 = qJD(1) * t118;
t192 = qJD(1) * t121;
t191 = qJD(1) * t144;
t190 = qJD(2) * t142;
t189 = qJD(2) * t145;
t188 = qJD(3) * t136;
t187 = qJD(3) * t142;
t186 = qJD(3) * t144;
t185 = qJD(5) * t105;
t184 = qJD(5) * t136;
t183 = t132 * qJD(4);
t182 = t141 * qJD(4);
t181 = t142 * qJD(1);
t180 = t144 * qJD(4);
t179 = t145 * qJD(1);
t174 = t115 * t181;
t173 = t117 * t181;
t172 = qJD(1) * t218;
t171 = t105 * t194;
t170 = t141 * t181;
t169 = t144 * t181;
t167 = t136 / 0.2e1 - t105 / 0.2e1;
t165 = pkin(4) * t177;
t122 = t178 * t142;
t163 = t178 * t144;
t123 = t178 * t145;
t162 = t120 / 0.2e1 - pkin(3) / 0.2e1 - t118 / 0.2e1;
t161 = t115 * t241;
t147 = t152 * t120;
t13 = t167 * t117 + t147;
t1 = t13 + t176;
t47 = t176 + t215;
t160 = -t1 * qJD(1) + t47 * qJD(3);
t158 = -t175 + t227;
t3 = t151 * t120 + t158;
t48 = t175 - t216;
t159 = t3 * qJD(1) + t48 * qJD(3);
t157 = -t189 + t193;
t156 = t190 + t192;
t155 = qJD(3) * t121 + t190;
t43 = t162 * t141;
t154 = qJD(1) * t43 + qJD(3) * t231;
t44 = t162 * t144;
t153 = pkin(3) * t186 + qJD(1) * t44;
t20 = -qJD(1) * t26 + t219;
t18 = t115 * t194 - t219;
t14 = -t167 * t115 + t226;
t150 = -t14 * qJD(1) - t115 * t188;
t149 = -t13 * qJD(1) + t117 * t188;
t16 = -t215 / 0.2e1 + t217 / 0.2e1 + t147;
t134 = t141 * t180;
t108 = (-t186 + t191) * t141;
t78 = t142 * t163 - t145 * t182;
t77 = t141 * t122 + t145 * t180;
t55 = -t208 / 0.2e1 - t148;
t54 = t209 / 0.2e1 + t146;
t46 = t168 + (pkin(3) + t118) * t232;
t45 = t231 / 0.2e1 + t118 * t141 / 0.2e1 + t212 / 0.2e1;
t37 = t51 * qJD(2);
t36 = t50 * qJD(2);
t22 = qJD(3) * t52 + t222;
t21 = qJD(3) * t49 + t223;
t15 = t226 + t227;
t8 = t177 * t26;
t6 = qJD(3) * t55 + t142 * t242 - t223;
t5 = qJD(3) * t54 + t142 * t241 - t222;
t4 = t158 + t226;
t2 = t16 - t176;
t7 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, t155, -t120 * qJD(3) + t189, t134, t183, 0, 0, 0, -t118 * t182 + t155 * t144, -t118 * t180 - t155 * t141, -t161, t243, 0, 0, 0, -qJD(4) * t29 - t155 * t115 - t117 * t185, -qJD(4) * t30 + t115 * t185 - t117 * t155; 0, 0, 0, 0, qJD(1), t195, 0, t181, t179, 0, 0, 0, 0, 0, t169, -t170, 0, 0, 0, 0, 0, t177 * t56 - t174, t177 * t53 - t173; 0, 0, 0, 0, 0, 0, 0, t164, -t178 * t120, -t134, -t183, 0, 0, 0, qJD(4) * t45 + t121 * t163, qJD(4) * t46 - t141 * t164, t8, -t243, 0, 0, 0, qJD(4) * t2 + qJD(5) * t16 - t115 * t164, qJD(4) * t4 + qJD(5) * t15 - t117 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t239, -t180, t182, 0, qJD(3) * t45 - t119 * t180 - t141 * t193, qJD(3) * t46 - t118 * t191 + t119 * t182, -t18, t11, t242, t241, 0, t2 * qJD(3) - t203 + t238, t4 * qJD(3) - t202 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t11, t242, t241, 0, t16 * qJD(3) - t171 + t238, t15 * qJD(3) + t172 + t237; 0, 0, 0, 0, -qJD(1), -t195, 0, -t122, -t123, 0, 0, 0, 0, 0, -t78, t77, 0, 0, 0, 0, 0, -qJD(3) * t57 - t177 * t50 + t174, -qJD(3) * t59 - t177 * t51 + t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t122, t123, 0, 0, 0, 0, 0, t78, -t77, 0, 0, 0, 0, 0, t115 * t187 + t177 * t55 - t221, t117 * t187 + t177 * t54 - t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t123 - t142 * t180, t144 * t123 + t142 * t182, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, -t156, t120 * qJD(1) - t189, -t134, -t183, 0, 0, 0, -qJD(4) * t43 - t156 * t144, -qJD(4) * t44 + t156 * t141, t8, -t243, 0, 0, 0, qJD(2) * t57 - qJD(4) * t1 - qJD(5) * t13 + t115 * t192, qJD(2) * t59 + qJD(4) * t3 - qJD(5) * t14 + t117 * t192; 0, 0, 0, 0, 0, 0, 0, -t181, -t179, 0, 0, 0, 0, 0, -t169, t170, 0, 0, 0, 0, 0, -t177 * t49 + t221, -t177 * t52 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t183, 0, 0, 0, -pkin(3) * t182, -pkin(3) * t180, -t161, t243, 0, 0, 0, qJD(4) * t47 + t117 * t184, qJD(4) * t48 - t115 * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t239, t180, -t182, 0, -pkin(7) * t180 - t154, pkin(7) * t182 - t153, -t20, -t11, -t242, -t241, 0, t160 + t236, t159 + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t11, -t242, -t241, 0, t149 + t236, t150 + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t239, 0, 0, 0, t43 * qJD(3) + t157 * t141, t44 * qJD(3) + t157 * t144, t18, -t11, 0, 0, 0, qJD(3) * t1 + t203 + t36, -qJD(3) * t3 + t202 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 * t179, -t144 * t179, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t239, 0, 0, 0, t154, t153, t20, t11, 0, 0, 0, -t160 + t197, -t159 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t225, -t143 * t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t165, -t143 * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t11, 0, 0, 0, qJD(3) * t13 + t171 + t36, qJD(3) * t14 - t172 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t11, 0, 0, 0, -t149 + t197, -t150 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t224, t143 * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
