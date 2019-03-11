% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:19
% EndTime: 2019-03-09 03:16:30
% DurationCPUTime: 3.43s
% Computational Cost: add. (6498->357), mult. (16986->466), div. (0->0), fcn. (13185->8), ass. (0->166)
t171 = sin(pkin(10));
t173 = cos(pkin(10));
t172 = sin(pkin(9));
t174 = cos(pkin(9));
t176 = sin(qJ(3));
t249 = cos(qJ(3));
t150 = t249 * t172 + t176 * t174;
t253 = t150 * qJD(1);
t118 = qJD(3) * t171 + t173 * t253;
t175 = sin(qJ(5));
t177 = cos(qJ(5));
t126 = t171 * t253;
t256 = qJD(3) * t173 - t126;
t192 = t177 * t256;
t72 = -t118 * t175 + t192;
t269 = t72 ^ 2;
t209 = t249 * t174;
t163 = qJD(1) * t209;
t223 = t172 * t176;
t208 = qJD(1) * t223;
t137 = -t163 + t208;
t132 = qJD(5) + t137;
t268 = t132 * t72;
t143 = t150 * qJD(3);
t128 = qJD(1) * t143;
t149 = t177 * t171 + t175 * t173;
t220 = t177 * t173;
t147 = t171 * t175 - t220;
t214 = qJD(5) * t177;
t215 = qJD(5) * t175;
t255 = -t171 * t215 + t173 * t214;
t264 = -t147 * t137 + t255;
t200 = t149 * t128 + t264 * t132;
t260 = t177 * t118 + t175 * t256;
t236 = t253 * t260;
t267 = t200 + t236;
t250 = t260 ^ 2;
t160 = qJD(3) * t163;
t127 = -qJD(3) * t208 + t160;
t228 = t127 * t171;
t40 = t175 * (qJD(5) * t118 + t228) - qJD(5) * t192 - t127 * t220;
t266 = t40 - t268;
t141 = t149 * qJD(5);
t232 = t149 * t137 + t141;
t199 = -t147 * t128 - t232 * t132;
t237 = t253 * t72;
t265 = t199 + t237;
t262 = t150 * qJD(2);
t261 = t253 * qJD(3);
t243 = pkin(8) + qJ(4);
t154 = t243 * t171;
t156 = t243 * t173;
t115 = -t154 * t175 + t156 * t177;
t247 = pkin(8) * t173;
t101 = pkin(3) * t253 + qJ(4) * t137;
t244 = pkin(7) + qJ(2);
t155 = t244 * t172;
t151 = qJD(1) * t155;
t157 = t244 * t174;
t152 = qJD(1) * t157;
t186 = t249 * t151 + t176 * t152;
t56 = t173 * t101 + t171 * t186;
t39 = pkin(4) * t253 + t137 * t247 + t56;
t227 = t137 * t171;
t57 = t171 * t101 - t173 * t186;
t48 = pkin(8) * t227 + t57;
t259 = -qJD(4) * t149 - qJD(5) * t115 + t175 * t48 - t177 * t39;
t193 = -t154 * t177 - t156 * t175;
t258 = qJD(4) * t147 - qJD(5) * t193 + t175 * t39 + t177 * t48;
t257 = t127 * t149;
t254 = -t249 * t155 - t176 * t157;
t252 = t147 * t40 - t232 * t260;
t185 = t209 - t223;
t168 = -pkin(2) * t174 - pkin(1);
t103 = -pkin(3) * t185 - qJ(4) * t150 + t168;
t116 = -t176 * t155 + t249 * t157;
t59 = t173 * t103 - t116 * t171;
t46 = -pkin(4) * t185 - t150 * t247 + t59;
t225 = t150 * t171;
t60 = t171 * t103 + t173 * t116;
t52 = -pkin(8) * t225 + t60;
t194 = t175 * t46 + t177 * t52;
t142 = t185 * qJD(3);
t79 = pkin(3) * t143 - qJ(4) * t142 - qJD(4) * t150;
t182 = t185 * qJD(2);
t86 = t254 * qJD(3) + t182;
t44 = -t171 * t86 + t173 * t79;
t25 = pkin(4) * t143 - t142 * t247 + t44;
t226 = t142 * t171;
t45 = t171 * t79 + t173 * t86;
t31 = -pkin(8) * t226 + t45;
t251 = -qJD(5) * t194 - t175 * t31 + t177 * t25;
t133 = t137 ^ 2;
t248 = pkin(5) * t128;
t98 = -qJD(3) * pkin(3) + qJD(4) + t186;
t62 = -pkin(4) * t256 + t98;
t16 = -pkin(5) * t72 - qJ(6) * t260 + t62;
t246 = t16 * t260;
t245 = t260 * t72;
t242 = qJ(6) * t253 + t258;
t241 = -pkin(5) * t253 + t259;
t107 = -t176 * t151 + t249 * t152;
t77 = -pkin(4) * t227 + t107;
t240 = -t232 * pkin(5) + t264 * qJ(6) + qJD(6) * t149 + t77;
t66 = pkin(3) * t128 - qJ(4) * t127 - qJD(4) * t253;
t73 = qJD(1) * t182 + (qJD(4) - t186) * qJD(3);
t33 = t171 * t66 + t173 * t73;
t100 = qJD(3) * qJ(4) + t107;
t153 = t168 * qJD(1) + qJD(2);
t83 = pkin(3) * t137 - qJ(4) * t253 + t153;
t51 = t173 * t100 + t171 * t83;
t50 = -t100 * t171 + t173 * t83;
t28 = pkin(4) * t137 - pkin(8) * t118 + t50;
t38 = pkin(8) * t256 + t51;
t9 = -t175 * t38 + t177 * t28;
t234 = qJD(6) - t9;
t231 = qJ(6) * t128;
t230 = t193 * t128;
t229 = t115 * t128;
t222 = t173 * t127;
t218 = t172 ^ 2 + t174 ^ 2;
t216 = qJD(3) * t176;
t213 = qJD(1) * qJD(2);
t167 = -t173 * pkin(4) - pkin(3);
t205 = qJD(3) * t249;
t32 = -t171 * t73 + t173 * t66;
t19 = pkin(4) * t128 - pkin(8) * t222 + t32;
t22 = -pkin(8) * t228 + t33;
t204 = t175 * t22 - t177 * t19 + t38 * t214 + t28 * t215;
t203 = t218 * qJD(1) ^ 2;
t76 = t262 * qJD(1) - t151 * t216 + t152 * t205;
t87 = -t155 * t216 + t157 * t205 + t262;
t41 = qJD(5) * t260 + t257;
t201 = -t149 * t41 + t264 * t72;
t198 = -t171 * t50 + t173 * t51;
t10 = t175 * t28 + t177 * t38;
t195 = -t175 * t52 + t177 * t46;
t55 = pkin(4) * t228 + t76;
t61 = pkin(4) * t226 + t87;
t190 = 0.2e1 * t218 * t213;
t88 = pkin(4) * t225 - t254;
t2 = t204 - t248;
t189 = t10 * t132 - t204;
t188 = t175 * t19 + t177 * t22 + t28 * t214 - t38 * t215;
t187 = t175 * t25 + t177 * t31 + t46 * t214 - t52 * t215;
t184 = -t127 * t254 + t98 * t142 + t76 * t150;
t181 = -pkin(3) * t127 - qJ(4) * t128 + (-qJD(4) + t98) * t137;
t5 = pkin(5) * t41 + qJ(6) * t40 - qJD(6) * t260 + t55;
t179 = t132 * t260 + t41;
t102 = pkin(5) * t147 - qJ(6) * t149 + t167;
t95 = t147 * t150;
t94 = t149 * t150;
t54 = t142 * t149 + t255 * t150;
t53 = t141 * t150 + t147 * t142;
t36 = pkin(5) * t260 - qJ(6) * t72;
t30 = t94 * pkin(5) + t95 * qJ(6) + t88;
t15 = -t40 - t268;
t14 = pkin(5) * t185 - t195;
t13 = -qJ(6) * t185 + t194;
t8 = qJ(6) * t132 + t10;
t7 = -pkin(5) * t132 + t234;
t6 = pkin(5) * t54 + qJ(6) * t53 + qJD(6) * t95 + t61;
t4 = -pkin(5) * t143 - t251;
t3 = qJ(6) * t143 - qJD(6) * t185 + t187;
t1 = qJD(6) * t132 + t188 + t231;
t11 = [0, 0, 0, 0, 0, t190, qJ(2) * t190, t127 * t150 + t142 * t253, t127 * t185 - t128 * t150 - t137 * t142 - t143 * t253, t142 * qJD(3), -t143 * qJD(3), 0, -qJD(3) * t87 + t128 * t168 + t143 * t153, -qJD(3) * t86 + t127 * t168 + t142 * t153, t59 * t128 + t44 * t137 + t50 * t143 + t184 * t171 - t185 * t32 - t256 * t87, t118 * t87 - t128 * t60 - t137 * t45 - t143 * t51 + t173 * t184 + t185 * t33, -t44 * t118 - t45 * t126 + (t45 * qJD(3) - t59 * t127 - t50 * t142 - t32 * t150) * t173 + (-t60 * t127 - t51 * t142 - t33 * t150) * t171, -t254 * t76 + t32 * t59 + t33 * t60 + t44 * t50 + t45 * t51 + t87 * t98, -t260 * t53 + t40 * t95, -t260 * t54 + t40 * t94 + t41 * t95 - t53 * t72, -t128 * t95 - t132 * t53 + t143 * t260 + t185 * t40, -t128 * t94 - t132 * t54 + t143 * t72 + t185 * t41, -t128 * t185 + t132 * t143, t195 * t128 + t251 * t132 + t9 * t143 + t185 * t204 + t88 * t41 + t62 * t54 + t55 * t94 - t61 * t72, -t10 * t143 - t194 * t128 - t187 * t132 + t185 * t188 + t260 * t61 - t88 * t40 - t62 * t53 - t55 * t95, -t128 * t14 - t132 * t4 - t143 * t7 + t16 * t54 + t185 * t2 + t30 * t41 + t5 * t94 - t6 * t72, -t1 * t94 - t13 * t41 - t14 * t40 - t2 * t95 + t260 * t4 + t3 * t72 - t53 * t7 - t54 * t8, -t1 * t185 + t128 * t13 + t132 * t3 + t143 * t8 + t16 * t53 - t260 * t6 + t30 * t40 + t5 * t95, t1 * t13 + t14 * t2 + t16 * t6 + t3 * t8 + t30 * t5 + t4 * t7; 0, 0, 0, 0, 0, -t203, -qJ(2) * t203, 0, 0, 0, 0, 0, 0.2e1 * t261, t160 + (-t137 - t208) * qJD(3), t173 * t128 - t171 * t133 + t253 * t256, -t118 * t253 - t128 * t171 - t133 * t173 (t171 * t118 + t173 * t256) * t137 + (-t171 ^ 2 - t173 ^ 2) * t127, t137 * t198 + t33 * t171 + t32 * t173 - t253 * t98, 0, 0, 0, 0, 0, t265, -t267, t265, t201 - t252, t267, t1 * t149 + t147 * t2 - t16 * t253 + t232 * t7 + t264 * t8; 0, 0, 0, 0, 0, 0, 0, t253 * t137, t253 ^ 2 - t133, t160 + (t137 - t208) * qJD(3), 0, 0, t107 * qJD(3) - t153 * t253 - t76, t153 * t137 - t185 * t213, t107 * t256 - t56 * t137 + t181 * t171 - t76 * t173 - t253 * t50, -t107 * t118 + t137 * t57 + t171 * t76 + t173 * t181 + t253 * t51, t56 * t118 + t57 * t126 + (-qJD(4) * t126 - t50 * t137 + t33 + (t173 * qJD(4) - t57) * qJD(3)) * t173 + (qJD(4) * t118 - t51 * t137 - t32) * t171, -pkin(3) * t76 - t107 * t98 - t50 * t56 - t51 * t57 + t198 * qJD(4) + (-t32 * t171 + t33 * t173) * qJ(4), -t149 * t40 + t260 * t264, t201 + t252, t200 - t236, t199 - t237, -t132 * t253, t259 * t132 + t55 * t147 + t167 * t41 + t232 * t62 - t253 * t9 + t72 * t77 + t230, t10 * t253 + t258 * t132 + t55 * t149 - t167 * t40 - t260 * t77 + t264 * t62 - t229, t102 * t41 + t241 * t132 + t147 * t5 + t232 * t16 + t240 * t72 + t253 * t7 + t230, -t1 * t147 - t115 * t41 + t149 * t2 + t193 * t40 - t232 * t8 - t241 * t260 - t242 * t72 + t264 * t7, t102 * t40 - t242 * t132 - t149 * t5 - t16 * t264 + t240 * t260 - t253 * t8 + t229, t1 * t115 + t102 * t5 - t240 * t16 - t193 * t2 - t241 * t7 - t242 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t137 + t228, t137 * t256 + t222, -t118 ^ 2 - t256 ^ 2, t118 * t50 - t256 * t51 + t76, 0, 0, 0, 0, 0, t179, -t266, t179, -t250 - t269, t266, -t260 * t7 - t72 * t8 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t250 - t269, t15, -t257 + (-qJD(5) + t132) * t260, t128, -t260 * t62 + t189, t132 * t9 - t62 * t72 - t188, t36 * t72 + t189 - t246 + 0.2e1 * t248, pkin(5) * t40 - qJ(6) * t41 + (-t10 + t8) * t260 - (t7 - t234) * t72, 0.2e1 * t231 + t16 * t72 + t36 * t260 + (0.2e1 * qJD(6) - t9) * t132 + t188, -pkin(5) * t2 + qJ(6) * t1 - t10 * t7 - t16 * t36 + t234 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245 - t261, t15, -t132 ^ 2 - t250, -t132 * t8 + t2 + t246;];
tauc_reg  = t11;
