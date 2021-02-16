% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:23
% EndTime: 2021-01-15 20:18:30
% DurationCPUTime: 2.40s
% Computational Cost: add. (3967->225), mult. (7586->284), div. (0->0), fcn. (8552->6), ass. (0->187)
t150 = qJD(2) + qJD(4);
t152 = sin(pkin(8));
t154 = sin(qJ(2));
t155 = cos(qJ(2));
t246 = cos(pkin(8));
t183 = t246 * t155;
t140 = -t152 * t154 + t183;
t232 = t152 * t155;
t141 = -t246 * t154 - t232;
t153 = sin(qJ(4));
t268 = cos(qJ(4));
t111 = t268 * t140 + t153 * t141;
t234 = t111 ^ 2;
t138 = t268 * t141;
t231 = t153 * t140;
t275 = -t138 + t231;
t283 = t275 ^ 2;
t284 = -t283 + t234;
t285 = t284 * qJD(1);
t280 = t111 * qJ(5);
t282 = t275 * pkin(4);
t58 = -t280 + t282;
t271 = -t111 / 0.2e1;
t272 = t275 / 0.2e1;
t147 = -pkin(2) * t155 - pkin(1);
t119 = -t140 * pkin(3) + t147;
t177 = -pkin(4) * t111 - qJ(5) * t275;
t47 = t119 + t177;
t281 = t111 * t47;
t256 = t275 * t47;
t279 = t150 * t111;
t278 = qJD(3) * t111;
t277 = t111 * qJD(1);
t208 = t111 * qJD(5);
t276 = t275 * qJD(1);
t179 = t246 * pkin(2) + pkin(3);
t267 = pkin(2) * t152;
t136 = t153 * t179 + t268 * t267;
t261 = -qJ(3) - pkin(6);
t184 = t152 * t261;
t181 = t154 * t184;
t118 = -t261 * t183 + t181;
t264 = t140 * pkin(7);
t162 = t118 + t264;
t158 = t268 * t162;
t178 = t261 * t246;
t171 = t154 * t178;
t116 = t261 * t232 + t171;
t262 = t141 * pkin(7);
t161 = t116 + t262;
t159 = t153 * t161;
t156 = t159 / 0.2e1 + t158 / 0.2e1;
t115 = t155 * t178 - t181;
t83 = t115 - t264;
t194 = t268 * t83;
t117 = t155 * t184 + t171;
t84 = t117 + t262;
t254 = t153 * t84;
t166 = -t254 / 0.2e1 + t194 / 0.2e1;
t29 = t156 + t166;
t228 = t29 * qJD(1);
t15 = qJD(2) * t136 + t228;
t135 = t153 * t267 - t268 * t179;
t160 = t153 * t162;
t80 = t268 * t161;
t157 = -t80 / 0.2e1 + t160 / 0.2e1;
t193 = t268 * t84;
t255 = t153 * t83;
t167 = -t255 / 0.2e1 - t193 / 0.2e1;
t26 = t157 - t167;
t247 = -t26 * qJD(1) - t135 * qJD(2);
t180 = -t138 / 0.2e1;
t106 = t180 + t138 / 0.2e1;
t28 = t156 - t166;
t52 = t158 + t159;
t197 = -t28 * qJD(2) + t106 * qJD(3) - t52 * qJD(4);
t274 = pkin(4) / 0.2e1;
t273 = -qJ(5) / 0.2e1;
t129 = qJ(5) + t136;
t270 = -t129 / 0.2e1;
t130 = -pkin(4) + t135;
t269 = -t130 / 0.2e1;
t263 = t141 * pkin(3);
t149 = t154 * pkin(2);
t71 = 0.2e1 * t180 + t231;
t260 = t29 * qJD(2) - t71 * qJD(3);
t49 = -t194 + t254;
t259 = -t49 * qJD(2) - t28 * qJD(4);
t121 = t149 - t263;
t48 = t121 + t58;
t9 = -t111 * t48 + t256;
t258 = qJD(1) * t9;
t257 = qJD(2) * pkin(2);
t50 = t160 - t80;
t51 = t193 + t255;
t3 = t47 * t48 + t49 * t50 + t51 * t52;
t253 = t3 * qJD(1);
t4 = t47 * t58;
t252 = t4 * qJD(1);
t6 = (t49 - t52) * t275 + (t50 + t51) * t111;
t251 = t6 * qJD(1);
t185 = t135 / 0.2e1 + t269;
t186 = t270 + t136 / 0.2e1;
t163 = -t111 * t185 + t186 * t275;
t169 = pkin(4) * t271 + t273 * t275;
t8 = t163 - t169;
t250 = t8 * qJD(1);
t249 = -qJD(3) * t275 - t29 * qJD(4);
t10 = -t275 * t48 - t281;
t245 = qJD(1) * t10;
t11 = t111 * t52 + t275 * t50;
t244 = qJD(1) * t11;
t12 = -t111 * t58 + t256;
t243 = qJD(1) * t12;
t13 = -t275 * t58 - t281;
t242 = qJD(1) * t13;
t38 = -t111 * t121 + t119 * t275;
t241 = qJD(1) * t38;
t39 = t111 * t119 + t121 * t275;
t240 = qJD(1) * t39;
t56 = t116 * t141 + t118 * t140;
t239 = qJD(1) * t56;
t85 = -t140 * t149 - t141 * t147;
t238 = qJD(1) * t85;
t86 = t140 * t147 - t141 * t149;
t237 = qJD(1) * t86;
t148 = t149 / 0.2e1;
t182 = t148 - t263 / 0.2e1;
t16 = (t273 + t270) * t111 + (t274 + t269) * t275 + t182;
t230 = t16 * qJD(1);
t30 = (t115 + t118) * t141 + (-t116 + t117) * t140;
t227 = t30 * qJD(1);
t32 = 0.2e1 * t272 * pkin(4) + 0.2e1 * t271 * qJ(5);
t226 = t32 * qJD(1);
t33 = t234 + t283;
t225 = t33 * qJD(1);
t35 = t115 * t116 + t117 * t118 + t147 * t149;
t223 = t35 * qJD(1);
t220 = t71 * qJD(1);
t165 = (t152 * t140 / 0.2e1 + t246 * t141 / 0.2e1) * pkin(2);
t82 = -t149 / 0.2e1 + t165;
t219 = t82 * qJD(1);
t217 = qJD(1) * t119;
t216 = qJD(1) * t155;
t212 = qJD(4) * t119;
t211 = t106 * qJD(1);
t91 = t106 * qJD(4);
t210 = t283 * qJD(1);
t114 = t140 ^ 2 + t141 ^ 2;
t206 = t114 * qJD(1);
t204 = t135 * qJD(4);
t203 = t140 * qJD(1);
t202 = t141 * qJD(1);
t145 = -t154 ^ 2 + t155 ^ 2;
t201 = t145 * qJD(1);
t200 = t154 * qJD(2);
t199 = t155 * qJD(2);
t198 = -t204 + qJD(5);
t196 = pkin(1) * t154 * qJD(1);
t195 = pkin(1) * t216;
t192 = t47 * t276;
t191 = t111 * t217;
t190 = t275 * t217;
t189 = t111 * t276;
t188 = t275 * t277;
t187 = t154 * t216;
t164 = -t185 * t52 + t186 * t50;
t170 = t51 * t273 + t49 * t274;
t2 = t164 + t170;
t57 = -t129 * t135 + t130 * t136;
t176 = t2 * qJD(1) + t57 * qJD(2);
t175 = qJD(2) * t26 + t278;
t27 = t157 + t167;
t174 = -qJD(2) * t27 - qJD(4) * t50;
t173 = qJD(2) * t51 - qJD(4) * t27;
t54 = qJD(2) * t275 + qJD(4) * t71;
t172 = -qJD(4) * t26 + t278;
t151 = qJ(5) * qJD(5);
t144 = t150 * qJ(5);
t126 = t136 * qJD(4);
t120 = t129 * qJD(5);
t81 = t148 + t165;
t17 = t111 * t129 / 0.2e1 + t130 * t272 + t282 / 0.2e1 - t280 / 0.2e1 + t182;
t14 = -t126 - t15;
t7 = t163 + t169;
t1 = t164 - t170;
t5 = [0, 0, 0, t154 * t199, t145 * qJD(2), 0, 0, 0, -pkin(1) * t200, -pkin(1) * t199, t85 * qJD(2), t86 * qJD(2), qJD(2) * t30 + qJD(3) * t114, qJD(2) * t35 + qJD(3) * t56, t279 * t275, t150 * t284, 0, 0, 0, qJD(2) * t38 + t212 * t275, qJD(2) * t39 + t111 * t212, qJD(2) * t9 + qJD(4) * t12 + t208 * t275, qJD(2) * t6 + qJD(3) * t33, qJD(2) * t10 + qJD(4) * t13 + qJD(5) * t283, qJD(2) * t3 + qJD(3) * t11 + qJD(4) * t4 - qJD(5) * t256; 0, 0, 0, t187, t201, t199, -t200, 0, -pkin(6) * t199 - t196, pkin(6) * t200 - t195, qJD(2) * t115 + t238, -qJD(2) * t117 + t237, t227 + (-t246 * t140 + t141 * t152) * t257, t223 + (t246 * t115 + t117 * t152) * t257 + t81 * qJD(3), t188, t285, t279, -t54, 0, t241 + t259, -t173 + t240, t258 + t259, t251 + (t111 * t130 - t129 * t275) * qJD(2) + t7 * qJD(4) + t208, t173 + t245, t253 + (t129 * t51 + t130 * t49) * qJD(2) + t17 * qJD(3) + t1 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, qJD(2) * t81 + t239, 0, 0, 0, 0, 0, t91, 0, t91, t225, 0, qJD(2) * t17 - qJD(5) * t106 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t285, t279, -qJD(2) * t71 - qJD(4) * t275, 0, t190 + t197, -t174 + t191, t197 + t243, t7 * qJD(2) + qJD(4) * t177 + t208, t174 + t242, t252 + t1 * qJD(2) + (-pkin(4) * t52 - qJ(5) * t50) * qJD(4) + t52 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t279, t210, -t192 - t197; 0, 0, 0, -t187, -t201, 0, 0, 0, t196, t195, qJD(3) * t141 - t238, -qJD(3) * t140 - t237, -t227, qJD(3) * t82 - t223, -t188, -t285, 0, -t91, 0, -t241 + t249, -t172 - t240, t249 - t258, qJD(4) * t8 - t251, t172 - t245, -qJD(3) * t16 + qJD(4) * t2 + qJD(5) * t29 - t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t204, -t126, 0, t198, qJD(4) * t57 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t203, 0, t219, 0, 0, 0, 0, 0, -t276, -t277, -t276, 0, t277, -t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, 0, t14, t204 - t247, t14, t250, t198 + t247, (-pkin(4) * t136 - qJ(5) * t135) * qJD(4) + t120 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t129 * t150 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 * qJD(2), t140 * qJD(2), -t206, -qJD(2) * t82 - t239, 0, 0, 0, 0, 0, t54, t279, t54, -t225, -t279, qJD(2) * t16 + qJD(4) * t32 - qJD(5) * t71 - t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t203, 0, -t219, 0, 0, 0, 0, 0, t276, t277, t276, 0, -t277, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t277, t220, 0, -t277, t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t285, 0, t106 * qJD(2), 0, -t190 + t260, -t175 - t191, -t243 + t260, -qJD(2) * t8, t175 - t242, -qJD(2) * t2 - qJD(3) * t32 - t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, t15, t247, t15, -t250, qJD(5) - t247, t151 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, -t277, -t220, 0, t277, -t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, -t210, t192 - t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -qJ(5) * qJD(4) - qJD(2) * t129 - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
