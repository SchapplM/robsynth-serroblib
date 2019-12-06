% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:58
% EndTime: 2019-12-05 16:28:08
% DurationCPUTime: 2.18s
% Computational Cost: add. (2041->201), mult. (5060->353), div. (0->0), fcn. (5789->10), ass. (0->181)
t152 = sin(pkin(10));
t277 = t152 / 0.2e1;
t154 = sin(qJ(5));
t158 = cos(qJ(3));
t265 = -qJ(4) - pkin(7);
t138 = t265 * t158;
t155 = sin(qJ(3));
t244 = t152 * t155;
t256 = cos(pkin(10));
t165 = -t256 * t138 + t265 * t244;
t250 = t165 * t154;
t157 = cos(qJ(5));
t249 = t165 * t157;
t174 = t256 * t158 - t244;
t127 = t174 ^ 2;
t195 = t256 * t155;
t243 = t152 * t158;
t131 = t195 + t243;
t128 = t131 ^ 2;
t276 = -t128 - t127;
t150 = t154 ^ 2;
t151 = t157 ^ 2;
t142 = t151 - t150;
t83 = t154 * t131;
t194 = 0.2e1 * t157 * t83;
t171 = qJD(2) * t194 - qJD(3) * t142;
t153 = sin(pkin(5));
t156 = sin(qJ(2));
t242 = t153 * t156;
t257 = cos(pkin(5));
t124 = t155 * t242 - t257 * t158;
t125 = t257 * t155 + t158 * t242;
t175 = -t152 * t124 + t256 * t125;
t275 = t175 / 0.2e1;
t271 = -t131 / 0.2e1;
t268 = -t157 / 0.2e1;
t267 = t157 / 0.2e1;
t266 = t155 * pkin(3);
t264 = qJD(3) * pkin(3);
t86 = pkin(4) * t131 - pkin(8) * t174 + t266;
t263 = t154 * t86;
t262 = t157 * t86;
t159 = cos(qJ(2));
t241 = t153 * t159;
t56 = t154 * t175 + t157 * t241;
t261 = t56 * t174;
t57 = -t154 * t241 + t157 * t175;
t260 = t57 * t174;
t69 = t256 * t124 + t125 * t152;
t259 = t69 * t131;
t258 = t175 * t174;
t217 = t128 - t127;
t58 = t217 * t154;
t255 = qJD(2) * t58;
t59 = t276 * t154;
t254 = qJD(2) * t59;
t60 = t217 * t157;
t253 = qJD(2) * t60;
t88 = t276 * t157;
t252 = qJD(2) * t88;
t10 = (-t175 / 0.2e1 + t275) * t131;
t251 = t10 * qJD(2);
t104 = t131 * t241;
t248 = t104 * t154;
t247 = t104 * t157;
t246 = t174 * t131;
t245 = t131 * t157;
t105 = t174 * t241;
t240 = t154 * t105;
t239 = t157 * t105;
t21 = -t153 ^ 2 * t156 * t159 + t69 * t104 + t105 * t175;
t237 = t21 * qJD(1);
t209 = -t256 / 0.2e1;
t168 = t131 * t209 + t174 * t277;
t67 = (-t155 / 0.2e1 + t168) * pkin(3);
t236 = t67 * qJD(2);
t80 = t154 * t174;
t235 = t80 * qJD(2);
t234 = t83 * qJD(2);
t85 = t157 * t174;
t233 = t85 * qJD(2);
t232 = t276 * qJD(2);
t230 = qJD(2) * t131;
t229 = qJD(2) * t156;
t228 = qJD(2) * t157;
t227 = qJD(2) * t158;
t226 = qJD(3) * t154;
t225 = qJD(3) * t157;
t224 = qJD(4) * t157;
t223 = qJD(5) * t154;
t222 = qJD(5) * t157;
t126 = t195 / 0.2e1 + t243 / 0.2e1;
t221 = t126 * qJD(2);
t143 = -t155 ^ 2 + t158 ^ 2;
t220 = t143 * qJD(2);
t219 = t155 * qJD(3);
t218 = t158 * qJD(3);
t216 = pkin(2) * t155 * qJD(2);
t215 = pkin(2) * t227;
t213 = t154 * t242;
t212 = t155 * t241;
t211 = -t259 / 0.2e1;
t210 = t69 * t267;
t148 = -pkin(3) * t158 - pkin(2);
t208 = t174 * t230;
t207 = t151 * t230;
t206 = t154 * t225;
t205 = t174 * t222;
t204 = qJD(3) * t246;
t203 = qJD(2) * t241;
t202 = t154 * t222;
t201 = t155 * t227;
t200 = t131 * t228;
t199 = t242 / 0.2e1;
t198 = -t240 / 0.2e1;
t197 = -t239 / 0.2e1;
t101 = -t138 * t152 - t265 * t195;
t193 = qJD(2) * t174 - qJD(5);
t191 = qJD(3) * t194;
t166 = t131 * t275;
t161 = t166 * t154 + t56 * t271;
t3 = t247 / 0.2e1 + t161;
t170 = -pkin(4) * t174 - pkin(8) * t131 + t148;
t32 = -t157 * t170 + t250;
t7 = (-t32 + t250) * t131 - t262 * t174;
t189 = t3 * qJD(1) + t7 * qJD(2);
t160 = t166 * t157 + t57 * t271;
t6 = -t248 / 0.2e1 + t160;
t33 = t154 * t170 + t249;
t8 = (-t33 + t249) * t131 + t263 * t174;
t188 = t6 * qJD(1) + t8 * qJD(2);
t169 = t104 * t209 + t105 * t277;
t1 = (t212 / 0.2e1 + t169) * pkin(3);
t24 = t148 * t266;
t187 = -t1 * qJD(1) + t24 * qJD(2);
t186 = t10 * qJD(1);
t180 = t199 + t211;
t12 = t198 - t260 / 0.2e1 + t180 * t157;
t23 = t101 * t245 + t174 * t33;
t185 = qJD(1) * t12 - qJD(2) * t23;
t13 = t197 + t261 / 0.2e1 - t180 * t154;
t22 = -t101 * t83 - t174 * t32;
t184 = -qJD(1) * t13 + qJD(2) * t22;
t25 = -t258 / 0.2e1 + t180;
t31 = t101 * t131 + t165 * t174;
t183 = -qJD(1) * t25 + qJD(2) * t31;
t146 = pkin(3) * t152 + pkin(8);
t147 = -t256 * pkin(3) - pkin(4);
t182 = -t131 * t146 + t147 * t174;
t181 = t193 * t157;
t179 = -t146 * t174 / 0.2e1 + t147 * t271;
t79 = (t150 / 0.2e1 - t151 / 0.2e1) * t131;
t178 = -qJD(2) * t79 + t206;
t177 = t131 * t181;
t176 = qJD(5) * t126 - t208;
t173 = t128 * t154 * t228 + qJD(3) * t79;
t87 = t142 * t128;
t172 = qJD(2) * t87 + t191;
t167 = t86 / 0.2e1 + t179;
t19 = t167 * t157;
t164 = qJD(2) * t19 - t147 * t226;
t17 = t167 * t154;
t163 = -qJD(2) * t17 - t147 * t225;
t123 = t126 * qJD(3);
t122 = t131 * t225;
t76 = t80 * qJD(5);
t75 = t79 * qJD(5);
t66 = t266 / 0.2e1 + t168 * pkin(3);
t65 = -t223 + t235;
t30 = -t268 * t69 + t210;
t29 = t69 * t154;
t26 = t258 / 0.2e1 + t259 / 0.2e1 + t199;
t20 = t262 / 0.2e1 - t179 * t157 + t101 * t154;
t18 = -t263 / 0.2e1 + t179 * t154 + (t267 - t268) * t101;
t15 = t260 / 0.2e1 + t131 * t210 + t198 + t157 * t199;
t14 = -t261 / 0.2e1 + t154 * t211 + t197 - t213 / 0.2e1;
t9 = t10 * qJD(3);
t5 = t248 / 0.2e1 + t160;
t4 = -t247 / 0.2e1 + t161;
t2 = (-t212 / 0.2e1 + t169) * pkin(3);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t153 * t229, -t203, 0, 0, 0, 0, 0, (-t156 * t227 - t159 * t219) * t153, (t155 * t229 - t159 * t218) * t153, (t104 * t131 + t105 * t174) * qJD(2) + t9, t237 + (t101 * t104 + t105 * t165 + t148 * t242) * qJD(2) + t2 * qJD(3) + t26 * qJD(4), 0, 0, 0, 0, 0, (-(t157 * t242 - t240) * t174 + t104 * t83) * qJD(2) + t4 * qJD(3) + t15 * qJD(5), ((t213 + t239) * t174 + t104 * t245) * qJD(2) + t5 * qJD(3) + t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 * qJD(3) - t155 * t203, t124 * qJD(3) - t158 * t203, t251, t2 * qJD(2) + (-t152 * t69 - t175 * t256) * t264, 0, 0, 0, 0, 0, qJD(2) * t4 + qJD(5) * t29 - t175 * t225, qJD(2) * t5 + qJD(5) * t30 + t175 * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t15 + qJD(3) * t29 - qJD(5) * t57, qJD(2) * t14 + qJD(3) * t30 + qJD(5) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -qJD(3) * t1 - qJD(4) * t25 - t237, 0, 0, 0, 0, 0, qJD(3) * t3 - qJD(5) * t12, qJD(3) * t6 - qJD(5) * t13; 0, 0, 0, 0, t155 * t218, t143 * qJD(3), 0, 0, 0, -pkin(2) * t219, -pkin(2) * t218, -qJD(4) * t276, qJD(3) * t24 + qJD(4) * t31, -t128 * t202 + t151 * t204, -qJD(5) * t87 - t174 * t191, qJD(3) * t60 + t223 * t246, -qJD(3) * t58 + t131 * t205, -t204, qJD(3) * t7 - qJD(4) * t59 + qJD(5) * t23, qJD(3) * t8 - qJD(4) * t88 + qJD(5) * t22; 0, 0, 0, 0, t201, t220, t218, -t219, 0, -pkin(7) * t218 - t216, pkin(7) * t219 - t215, (-t131 * t152 - t174 * t256) * t264 + t186, (-t101 * t152 - t165 * t256) * t264 + t66 * qJD(4) + t187, -t75 - (-t206 - t207) * t174, -0.2e1 * t131 * t202 - t171 * t174, t131 * t226 + t253, t122 - t255, t176, (t154 * t182 - t249) * qJD(3) + t20 * qJD(5) + t189, (t157 * t182 + t250) * qJD(3) + t18 * qJD(5) + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, qJD(3) * t66 + t183, 0, 0, 0, 0, 0, -t254, -t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, -t172, t193 * t83, t177, t123, qJD(3) * t20 - qJD(5) * t33 - t185, qJD(3) * t18 + qJD(5) * t32 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, qJD(2) * t1, 0, 0, 0, 0, 0, -qJD(2) * t3, -qJD(2) * t6; 0, 0, 0, 0, -t201, -t220, 0, 0, 0, t216, t215, -t186, qJD(4) * t67 - t187, -t174 * t207 - t75, 0.2e1 * t154 * t177, -qJD(5) * t85 - t253, t76 + t255, -t176, -qJD(5) * t19 - t131 * t224 - t189, qJD(4) * t83 + qJD(5) * t17 - t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t142 * qJD(5), 0, 0, 0, t147 * t223, t147 * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, 0, 0, 0, 0, 0, -t200, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, -t171, t222 - t233, t65, -t221, -t146 * t222 - t164, t146 * t223 - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, -qJD(3) * t67 - t183, 0, 0, 0, 0, 0, t122 + t76 + t254, -qJD(3) * t83 + t205 + t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, 0, 0, 0, 0, 0, t200, -t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t12, qJD(2) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t172, qJD(3) * t85 - t154 * t208, -qJD(3) * t80 - t174 * t200, t123, qJD(3) * t19 - qJD(4) * t80 + t185, -qJD(3) * t17 - t174 * t224 - t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t171, t233, -t235, t221, t164, t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, -t174 * t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
