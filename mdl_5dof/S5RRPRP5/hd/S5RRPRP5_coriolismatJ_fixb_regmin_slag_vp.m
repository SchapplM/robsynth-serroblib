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
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 19:54:59
% EndTime: 2019-12-31 19:55:05
% DurationCPUTime: 2.13s
% Computational Cost: add. (3925->217), mult. (7494->268), div. (0->0), fcn. (8448->6), ass. (0->182)
t148 = qJD(2) + qJD(4);
t150 = sin(pkin(8));
t152 = sin(qJ(2));
t153 = cos(qJ(2));
t241 = cos(pkin(8));
t181 = t241 * t153;
t138 = -t150 * t152 + t181;
t232 = t150 * t153;
t139 = -t241 * t152 - t232;
t151 = sin(qJ(4));
t263 = cos(qJ(4));
t109 = t263 * t138 + t151 * t139;
t236 = t109 ^ 2;
t136 = t263 * t139;
t231 = t151 * t138;
t270 = -t136 + t231;
t278 = t270 ^ 2;
t279 = -t278 + t236;
t280 = t279 * qJD(1);
t275 = t109 * qJ(5);
t277 = t270 * pkin(4);
t58 = -t275 + t277;
t266 = -t109 / 0.2e1;
t267 = t270 / 0.2e1;
t190 = -t153 * pkin(2) - pkin(1);
t117 = -t138 * pkin(3) + t190;
t175 = -pkin(4) * t109 - qJ(5) * t270;
t47 = t117 + t175;
t276 = t47 * t109;
t247 = t47 * t270;
t209 = t109 * qJD(4);
t274 = t109 * qJD(2) + t209;
t273 = t109 * qJD(1);
t272 = t109 * qJD(3);
t208 = t109 * qJD(5);
t271 = t270 * qJD(1);
t177 = t241 * pkin(2) + pkin(3);
t262 = pkin(2) * t150;
t134 = t151 * t177 + t263 * t262;
t256 = -qJ(3) - pkin(6);
t182 = t150 * t256;
t179 = t152 * t182;
t116 = -t256 * t181 + t179;
t259 = t138 * pkin(7);
t160 = t116 + t259;
t156 = t263 * t160;
t176 = t256 * t241;
t169 = t152 * t176;
t114 = t256 * t232 + t169;
t257 = t139 * pkin(7);
t159 = t114 + t257;
t157 = t151 * t159;
t154 = t157 / 0.2e1 + t156 / 0.2e1;
t113 = t153 * t176 - t179;
t83 = t113 - t259;
t193 = t263 * t83;
t115 = t153 * t182 + t169;
t84 = t115 + t257;
t250 = t151 * t84;
t164 = -t250 / 0.2e1 + t193 / 0.2e1;
t29 = t154 + t164;
t228 = t29 * qJD(1);
t15 = t134 * qJD(2) + t228;
t133 = t151 * t262 - t263 * t177;
t158 = t151 * t160;
t80 = t263 * t159;
t155 = -t80 / 0.2e1 + t158 / 0.2e1;
t192 = t263 * t84;
t251 = t151 * t83;
t165 = -t251 / 0.2e1 - t192 / 0.2e1;
t26 = t155 - t165;
t242 = -t26 * qJD(1) - t133 * qJD(2);
t178 = -t136 / 0.2e1;
t104 = t178 + t136 / 0.2e1;
t28 = t154 - t164;
t52 = t156 + t157;
t196 = -t28 * qJD(2) + t104 * qJD(3) - t52 * qJD(4);
t269 = pkin(4) / 0.2e1;
t268 = -qJ(5) / 0.2e1;
t127 = qJ(5) + t134;
t265 = -t127 / 0.2e1;
t128 = -pkin(4) + t133;
t264 = -t128 / 0.2e1;
t258 = t139 * pkin(3);
t147 = t152 * pkin(2);
t71 = 0.2e1 * t178 + t231;
t255 = t29 * qJD(2) - t71 * qJD(3);
t254 = -qJD(3) * t270 - t29 * qJD(4);
t49 = -t193 + t250;
t253 = -t49 * qJD(2) - t28 * qJD(4);
t252 = qJD(2) * pkin(2);
t119 = t147 - t258;
t48 = t119 + t58;
t50 = t158 - t80;
t51 = t192 + t251;
t3 = t47 * t48 + t50 * t49 + t52 * t51;
t249 = t3 * qJD(1);
t4 = t47 * t58;
t248 = t4 * qJD(1);
t6 = (t49 - t52) * t270 + (t50 + t51) * t109;
t246 = t6 * qJD(1);
t183 = t133 / 0.2e1 + t264;
t184 = t265 + t134 / 0.2e1;
t161 = -t109 * t183 + t184 * t270;
t167 = pkin(4) * t266 + t268 * t270;
t8 = t161 - t167;
t245 = t8 * qJD(1);
t9 = -t109 * t48 + t247;
t244 = t9 * qJD(1);
t10 = -t270 * t48 - t276;
t240 = t10 * qJD(1);
t11 = t52 * t109 + t270 * t50;
t235 = t11 * qJD(1);
t12 = -t109 * t58 + t247;
t234 = t12 * qJD(1);
t13 = -t270 * t58 - t276;
t233 = t13 * qJD(1);
t146 = t147 / 0.2e1;
t180 = t146 - t258 / 0.2e1;
t16 = (t268 + t265) * t109 + (t269 + t264) * t270 + t180;
t230 = t16 * qJD(1);
t30 = (t113 + t116) * t139 + (-t114 + t115) * t138;
t227 = t30 * qJD(1);
t32 = 0.2e1 * t267 * pkin(4) + 0.2e1 * t266 * qJ(5);
t226 = t32 * qJD(1);
t33 = t236 + t278;
t225 = t33 * qJD(1);
t35 = t114 * t113 + t116 * t115 + t190 * t147;
t223 = t35 * qJD(1);
t38 = -t109 * t119 + t117 * t270;
t222 = t38 * qJD(1);
t39 = t117 * t109 + t119 * t270;
t221 = t39 * qJD(1);
t56 = t114 * t139 + t116 * t138;
t219 = t56 * qJD(1);
t217 = t71 * qJD(1);
t163 = (t150 * t138 / 0.2e1 + t241 * t139 / 0.2e1) * pkin(2);
t82 = -t147 / 0.2e1 + t163;
t216 = t82 * qJD(1);
t214 = qJD(1) * t117;
t213 = qJD(1) * t153;
t212 = t104 * qJD(1);
t89 = t104 * qJD(4);
t211 = t278 * qJD(1);
t205 = t270 * qJD(4);
t112 = t138 ^ 2 + t139 ^ 2;
t204 = t112 * qJD(1);
t202 = t133 * qJD(4);
t144 = -t152 ^ 2 + t153 ^ 2;
t200 = t144 * qJD(1);
t199 = t152 * qJD(2);
t198 = t153 * qJD(2);
t197 = -t202 + qJD(5);
t195 = pkin(1) * t152 * qJD(1);
t194 = pkin(1) * t213;
t191 = t47 * t271;
t189 = t109 * t271;
t188 = t270 * t273;
t187 = t109 * t214;
t186 = t270 * t214;
t185 = t152 * t213;
t162 = -t183 * t52 + t184 * t50;
t168 = t51 * t268 + t49 * t269;
t2 = t162 + t168;
t57 = -t127 * t133 + t128 * t134;
t174 = t2 * qJD(1) + t57 * qJD(2);
t173 = t26 * qJD(2) + t272;
t27 = t155 + t165;
t172 = -t27 * qJD(2) - t50 * qJD(4);
t171 = t51 * qJD(2) - t27 * qJD(4);
t54 = qJD(2) * t270 + t71 * qJD(4);
t170 = -t26 * qJD(4) + t272;
t149 = qJ(5) * qJD(5);
t143 = t148 * qJ(5);
t124 = t134 * qJD(4);
t118 = t127 * qJD(5);
t81 = t146 + t163;
t17 = t109 * t127 / 0.2e1 + t128 * t267 + t277 / 0.2e1 - t275 / 0.2e1 + t180;
t14 = -t124 - t15;
t7 = t161 + t167;
t1 = t162 - t168;
t5 = [0, 0, 0, t152 * t198, t144 * qJD(2), 0, 0, 0, -pkin(1) * t199, -pkin(1) * t198, t30 * qJD(2) + t112 * qJD(3), t35 * qJD(2) + t56 * qJD(3), t274 * t270, t148 * t279, 0, 0, 0, t38 * qJD(2) + t117 * t205, t39 * qJD(2) + t117 * t209, t9 * qJD(2) + t12 * qJD(4) + t208 * t270, t6 * qJD(2) + t33 * qJD(3), t10 * qJD(2) + t13 * qJD(4) + qJD(5) * t278, t3 * qJD(2) + t11 * qJD(3) + t4 * qJD(4) - qJD(5) * t247; 0, 0, 0, t185, t200, t198, -t199, 0, -pkin(6) * t198 - t195, pkin(6) * t199 - t194, t227 + (-t241 * t138 + t139 * t150) * t252, t223 + (t241 * t113 + t115 * t150) * t252 + t81 * qJD(3), t188, t280, t274, -t54, 0, t222 + t253, -t171 + t221, t244 + t253, t246 + (t128 * t109 - t127 * t270) * qJD(2) + t7 * qJD(4) + t208, t171 + t240, t249 + (t51 * t127 + t49 * t128) * qJD(2) + t17 * qJD(3) + t1 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, t81 * qJD(2) + t219, 0, 0, 0, 0, 0, t89, 0, t89, t225, 0, t17 * qJD(2) - t104 * qJD(5) + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t280, t274, -t71 * qJD(2) - t205, 0, t186 + t196, -t172 + t187, t196 + t234, t7 * qJD(2) + qJD(4) * t175 + t208, t172 + t233, t248 + t1 * qJD(2) + (-t52 * pkin(4) - t50 * qJ(5)) * qJD(4) + t52 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t274, t211, -t191 - t196; 0, 0, 0, -t185, -t200, 0, 0, 0, t195, t194, -t227, t82 * qJD(3) - t223, -t188, -t280, 0, -t89, 0, -t222 + t254, -t170 - t221, -t244 + t254, t8 * qJD(4) - t246, t170 - t240, -t16 * qJD(3) + t2 * qJD(4) + t29 * qJD(5) - t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t202, -t124, 0, t197, t57 * qJD(4) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, 0, 0, 0, 0, 0, -t271, -t273, -t271, 0, t273, -t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, 0, t14, t202 - t242, t14, t245, t197 + t242, (-t134 * pkin(4) - t133 * qJ(5)) * qJD(4) + t118 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t127 * t148 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, -t82 * qJD(2) - t219, 0, 0, 0, 0, 0, t54, t274, t54, -t225, -t274, t16 * qJD(2) + t32 * qJD(4) - t71 * qJD(5) - t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216, 0, 0, 0, 0, 0, t271, t273, t271, 0, -t273, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, t273, t217, 0, -t273, t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t280, 0, t104 * qJD(2), 0, -t186 + t255, -t173 - t187, -t234 + t255, -t8 * qJD(2), t173 - t233, -t2 * qJD(2) - t32 * qJD(3) - t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, 0, t15, t242, t15, -t245, qJD(5) - t242, t149 - t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, -t273, -t217, 0, t273, -t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, 0, -t211, t191 - t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -qJ(5) * qJD(4) - t127 * qJD(2) - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
