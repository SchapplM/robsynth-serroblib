% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:06
% EndTime: 2019-12-31 18:52:12
% DurationCPUTime: 1.85s
% Computational Cost: add. (3376->221), mult. (6721->322), div. (0->0), fcn. (7177->6), ass. (0->198)
t149 = sin(pkin(8));
t259 = cos(qJ(3));
t189 = t259 * t149;
t150 = cos(pkin(8));
t152 = sin(qJ(3));
t226 = t152 * t150;
t131 = t189 + t226;
t151 = sin(qJ(4));
t153 = cos(qJ(4));
t228 = t151 * t153;
t180 = 0.2e1 * t131 * t228;
t129 = t152 * t149 - t259 * t150;
t257 = t129 * pkin(4);
t234 = t131 * t153;
t253 = pkin(6) + qJ(2);
t134 = t253 * t150;
t190 = t259 * t134;
t133 = t253 * t149;
t227 = t152 * t133;
t88 = t190 - t227;
t249 = t151 * t88;
t143 = -t150 * pkin(2) - pkin(1);
t175 = t129 * pkin(3) - t131 * pkin(7);
t64 = t143 + t175;
t46 = -t153 * t64 + t249;
t42 = -qJ(5) * t234 - t46;
t29 = t42 + t257;
t266 = -t29 + t42;
t252 = -qJ(5) - pkin(7);
t136 = t252 * t153;
t261 = t42 / 0.2e1;
t193 = t261 - t29 / 0.2e1;
t265 = t193 * t136;
t126 = t129 ^ 2;
t127 = t131 ^ 2;
t264 = -t127 - t126;
t197 = t127 - t126;
t123 = t189 / 0.2e1 + t226 / 0.2e1;
t87 = t259 * t133 + t152 * t134;
t237 = t87 * t151;
t71 = t153 * t129;
t255 = t131 * pkin(3);
t256 = t129 * pkin(7);
t85 = t255 + t256;
t79 = t153 * t85;
t35 = t131 * pkin(4) + qJ(5) * t71 + t237 + t79;
t262 = t35 / 0.2e1;
t148 = t153 ^ 2;
t260 = -t148 / 0.2e1;
t258 = pkin(4) * t151;
t254 = t153 * pkin(4);
t248 = t153 * t88;
t47 = t151 * t64 + t248;
t68 = t151 * t131;
t43 = -qJ(5) * t68 + t47;
t243 = t43 * t153;
t247 = t29 * t151;
t251 = -t247 / 0.2e1 + t243 / 0.2e1;
t250 = pkin(4) * qJD(4);
t66 = t151 * t129;
t78 = t151 * t85;
t84 = t87 * t153;
t45 = qJ(5) * t66 + t78 - t84;
t55 = pkin(4) * t68 + t87;
t56 = -pkin(4) * t66 + t88;
t3 = t29 * t35 + t43 * t45 + t55 * t56;
t246 = t3 * qJD(1);
t245 = t35 * t153;
t194 = t257 / 0.2e1;
t168 = t194 - t193;
t4 = t168 * t153;
t244 = t4 * qJD(1);
t242 = t45 * t151;
t241 = t55 * t151;
t173 = t43 * t151 + t29 * t153;
t6 = (t242 + t245) * t131 - t173 * t129;
t240 = t6 * qJD(1);
t7 = t266 * t68;
t239 = t7 * qJD(1);
t196 = pkin(4) * t234;
t8 = t55 * t196 + t266 * t43;
t238 = t8 * qJD(1);
t9 = t168 * t151;
t236 = t9 * qJD(1);
t13 = t55 * t131 + (-t243 + t247) * t129;
t235 = t13 * qJD(1);
t135 = t252 * t151;
t233 = t135 * t151;
t232 = t135 * t153;
t231 = t136 * t151;
t230 = t136 * t153;
t144 = -pkin(3) - t254;
t181 = t131 * t144 / 0.2e1;
t154 = (t230 / 0.2e1 + t233 / 0.2e1) * t129 + t181;
t167 = -t245 / 0.2e1 - t242 / 0.2e1;
t14 = t154 + t167;
t229 = t14 * qJD(1);
t16 = (-t46 + t249) * t131 + t79 * t129;
t225 = t16 * qJD(1);
t17 = (-t47 + t248) * t131 - t78 * t129;
t224 = t17 * qJD(1);
t18 = t173 * t131;
t223 = t18 * qJD(1);
t23 = t46 * t129 - t87 * t68;
t222 = t23 * qJD(1);
t24 = -t47 * t129 + t87 * t234;
t221 = t24 * qJD(1);
t147 = t151 ^ 2;
t176 = (-t147 / 0.2e1 + t260) * t131;
t51 = t176 - t123;
t220 = t51 * qJD(1);
t52 = t197 * t151;
t219 = t52 * qJD(1);
t53 = t264 * t151;
t218 = t53 * qJD(1);
t54 = t197 * t153;
t217 = t54 * qJD(1);
t216 = t197 * qJD(1);
t215 = t66 * qJD(1);
t214 = t68 * qJD(1);
t213 = t71 * qJD(1);
t118 = t147 * t129;
t119 = t148 * t129;
t74 = t118 + t119;
t212 = t74 * qJD(1);
t140 = t147 + t148;
t75 = t140 * t127;
t211 = t75 * qJD(1);
t77 = t264 * t153;
t210 = t77 * qJD(1);
t137 = t149 ^ 2 + t150 ^ 2;
t141 = t148 - t147;
t209 = qJD(2) * t153;
t208 = qJD(3) * t151;
t207 = qJD(3) * t153;
t206 = qJD(4) * t151;
t205 = qJD(4) * t153;
t204 = t123 * qJD(1);
t203 = t129 * qJD(1);
t120 = t129 * qJD(3);
t202 = t131 * qJD(1);
t201 = t131 * qJD(3);
t132 = t137 * qJ(2);
t200 = t132 * qJD(1);
t199 = t137 * qJD(1);
t198 = t140 * qJD(3);
t195 = pkin(4) * t206;
t192 = t78 / 0.2e1 - t84 / 0.2e1;
t188 = t131 * t206;
t187 = t131 * t205;
t186 = t129 * t202;
t185 = t129 * t201;
t184 = t151 * t205;
t183 = t151 * t207;
t182 = t153 * t202;
t179 = t151 * t194;
t178 = qJD(1) * t143 + qJD(2);
t177 = -qJD(4) - t203;
t174 = qJD(3) * t180;
t1 = t265 + (t262 - t144 * t234 / 0.2e1 - t241 / 0.2e1) * pkin(4);
t49 = t144 * t258;
t172 = -t1 * qJD(1) + t49 * qJD(3);
t155 = t227 / 0.2e1 - t190 / 0.2e1 + t179;
t157 = (t231 / 0.2e1 - t232 / 0.2e1) * t131 + t251;
t12 = t155 + t157;
t97 = -t230 - t233;
t171 = t12 * qJD(1) + t97 * qJD(3);
t170 = t177 * t153;
t169 = t256 / 0.2e1 + t255 / 0.2e1;
t161 = t169 * t153;
t21 = -t79 / 0.2e1 - t161;
t166 = pkin(3) * t208 - t21 * qJD(1);
t156 = t169 * t151 + t84 / 0.2e1;
t19 = t156 + t192;
t165 = pkin(3) * t207 - t19 * qJD(1);
t65 = (t147 / 0.2e1 + t260) * t131;
t164 = -t65 * qJD(1) + t183;
t163 = t131 * t170;
t162 = t123 * qJD(4) + t186;
t160 = t127 * qJD(1) * t228 + t65 * qJD(3);
t76 = t141 * t127;
t159 = t76 * qJD(1) + t174;
t158 = qJD(1) * t180 - t141 * qJD(3);
t117 = t123 * qJD(3);
t115 = t153 * t201;
t98 = (t182 + t208) * pkin(4);
t61 = t66 * qJD(4);
t60 = t65 * qJD(4);
t57 = -t206 - t215;
t50 = t176 + t123;
t22 = t237 + t79 / 0.2e1 - t161;
t20 = t156 - t192;
t15 = t154 - t167;
t11 = -t155 + t157;
t10 = t151 * t261 - t243 / 0.2e1 + t179 + t251;
t5 = (t193 + t194) * t153;
t2 = t181 * t254 - t265 + (t241 / 0.2e1 + t262) * pkin(4);
t25 = [0, 0, 0, 0, 0, t137 * qJD(2), t132 * qJD(2), -t185, -t197 * qJD(3), 0, 0, 0, t143 * t201, -t143 * t120, -t127 * t184 - t148 * t185, -t76 * qJD(4) + t129 * t174, t54 * qJD(3) - t129 * t188, -t52 * qJD(3) - t129 * t187, t185, -t53 * qJD(2) + t16 * qJD(3) + t24 * qJD(4), -t77 * qJD(2) + t17 * qJD(3) + t23 * qJD(4), -t6 * qJD(3) - t7 * qJD(4) + t75 * qJD(5), t13 * qJD(2) + t3 * qJD(3) + t8 * qJD(4) - t18 * qJD(5); 0, 0, 0, 0, 0, t199, t200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, -t210, 0, t15 * qJD(3) + t10 * qJD(4) + t50 * qJD(5) + t235; 0, 0, 0, 0, 0, 0, 0, -t186, -t216, -t120, -t201, 0, -t88 * qJD(3) + t143 * t202, t87 * qJD(3) - t143 * t203, -t60 + (-t148 * t202 - t183) * t129, (t118 - t119) * qJD(3) + (-qJD(4) + t203) * t180, t151 * t201 + t217, t115 - t219, t162, t225 + (t151 * t175 - t248) * qJD(3) + t22 * qJD(4), t224 + (t153 * t175 + t249) * qJD(3) + t20 * qJD(4), -t240 + (-t35 * t151 + t45 * t153 + (-t231 + t232) * t129) * qJD(3) + t5 * qJD(4), t246 + t15 * qJD(2) + (t35 * t135 - t45 * t136 + t56 * t144) * qJD(3) + t2 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, -t159, t177 * t68, t163, t117, t22 * qJD(3) - t47 * qJD(4) + t221, t20 * qJD(3) + t46 * qJD(4) + t222, pkin(4) * t188 + t5 * qJD(3) - t239, t10 * qJD(2) + t2 * qJD(3) - t250 * t43 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t50 * qJD(2) + t11 * qJD(3) - t223; 0, 0, 0, 0, 0, -t199, -t200, 0, 0, 0, 0, 0, t201, -t120, 0, 0, 0, 0, 0, t115 - t61 + t218, -t68 * qJD(3) - t129 * t205 + t210, t74 * qJD(3), -t14 * qJD(3) - t9 * qJD(4) + t51 * qJD(5) - t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t203, 0, 0, 0, 0, 0, t182, -t214, t212, -t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t170, 0, -t195 - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220; 0, 0, 0, 0, 0, 0, 0, t186, t216, 0, 0, 0, -t178 * t131, t178 * t129, t148 * t186 - t60, 0.2e1 * t151 * t163, t71 * qJD(4) - t217, -t61 + t219, -t162, t21 * qJD(4) - t131 * t209 - t225, t68 * qJD(2) + t19 * qJD(4) - t224, -t74 * qJD(2) - t4 * qJD(4) + t240, t14 * qJD(2) - t1 * qJD(4) + t12 * qJD(5) - t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t203, 0, 0, 0, 0, 0, -t182, t214, -t212, t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t141 * qJD(4), 0, 0, 0, -pkin(3) * t206, -pkin(3) * t205, t140 * qJD(5), t49 * qJD(4) + t97 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t158, t205 + t213, t57, -t204, -pkin(7) * t205 - t166, pkin(7) * t206 - t165, -pkin(4) * t205 - t244, t136 * t250 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t159, -t71 * qJD(3) + t151 * t186, t66 * qJD(3) + t129 * t182, t117, t66 * qJD(2) - t21 * qJD(3) - t221, -t19 * qJD(3) + t129 * t209 - t222, t4 * qJD(3) + t239, t9 * qJD(2) + t1 * qJD(3) - qJD(5) * t196 - t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, t153 * t203, 0, t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t158, -t213, t215, t204, t166, t165, t244, -qJD(5) * t258 - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, pkin(4) * t187 - t51 * qJD(2) - t12 * qJD(3) + t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t171 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t25;
