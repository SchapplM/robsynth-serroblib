% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:21
% EndTime: 2019-12-05 16:52:26
% DurationCPUTime: 1.64s
% Computational Cost: add. (1914->184), mult. (4468->259), div. (0->0), fcn. (4527->6), ass. (0->160)
t155 = qJD(3) + qJD(4);
t161 = cos(qJ(2));
t160 = cos(qJ(3));
t250 = cos(qJ(4));
t189 = t250 * t160;
t157 = sin(qJ(4));
t158 = sin(qJ(3));
t228 = t157 * t158;
t168 = t189 - t228;
t270 = t168 * t161;
t190 = t250 * t158;
t227 = t157 * t160;
t135 = t190 + t227;
t166 = -t227 / 0.2e1 - t190 / 0.2e1;
t62 = (t135 / 0.2e1 + t166) * t161;
t215 = t62 * qJD(1);
t258 = -pkin(7) - pkin(6);
t144 = t258 * t160;
t164 = -t250 * t144 + t258 * t228;
t267 = t155 * t164;
t269 = -t215 - t267;
t180 = -t157 * t144 - t258 * t190;
t268 = t155 * t180;
t196 = t250 / 0.2e1;
t159 = sin(qJ(2));
t117 = t135 * t159;
t266 = t155 * t117;
t119 = t168 * t159;
t265 = t155 * t119;
t152 = -t160 * pkin(3) - pkin(2);
t208 = qJD(2) * t152;
t183 = t270 / 0.2e1;
t226 = t157 * t161;
t184 = t226 / 0.2e1;
t251 = -t161 / 0.2e1;
t212 = t158 * t184 + t189 * t251;
t64 = t183 + t212;
t213 = t64 * qJD(1);
t264 = -t168 * t208 + t213;
t246 = t158 * pkin(3);
t60 = t135 * t246 + t152 * t168;
t263 = -t60 * qJD(2) + t213;
t182 = -t270 / 0.2e1;
t179 = t161 * t196;
t211 = -t158 * t226 / 0.2e1 + t160 * t179;
t63 = t182 + t211;
t214 = t63 * qJD(1);
t231 = t135 * qJ(5);
t248 = t168 * pkin(4);
t178 = -t231 - t248;
t70 = t152 + t178;
t241 = t70 * t168;
t83 = t135 * pkin(4) - qJ(5) * t168;
t26 = -t83 * t135 - t241;
t262 = -t26 * qJD(2) + t214;
t73 = t83 + t246;
t23 = -t73 * t135 - t241;
t261 = -t23 * qJD(2) + t214;
t259 = t164 * qJD(5);
t132 = t135 ^ 2;
t257 = -qJ(5) / 0.2e1;
t256 = t164 / 0.2e1;
t118 = t161 * t135;
t181 = -t118 / 0.2e1;
t255 = t118 / 0.2e1;
t247 = t157 * pkin(3);
t149 = qJ(5) + t247;
t254 = -t149 / 0.2e1;
t195 = t250 * pkin(3);
t151 = -t195 - pkin(4);
t253 = t151 / 0.2e1;
t252 = t157 / 0.2e1;
t249 = t119 * pkin(4);
t242 = pkin(3) * qJD(4);
t240 = t70 * t135;
t239 = t119 * qJD(5);
t234 = t117 * qJ(5);
t233 = t117 * t149;
t232 = t119 * t151;
t230 = t149 * t135;
t229 = t151 * t168;
t167 = t195 / 0.2e1 + t253 + pkin(4) / 0.2e1;
t171 = t254 + t247 / 0.2e1 + qJ(5) / 0.2e1;
t21 = t171 * t135 + t167 * t168;
t224 = t21 * qJD(2);
t22 = -t168 * t73 + t240;
t223 = t22 * qJD(2);
t25 = -t168 * t83 + t240;
t221 = t25 * qJD(2);
t27 = t117 * t118 + t119 * t270 - t161 * t159;
t219 = t27 * qJD(1);
t46 = t168 ^ 2 - t132;
t218 = t46 * qJD(2);
t59 = t152 * t135 - t168 * t246;
t217 = t59 * qJD(2);
t210 = t158 * t179 + t160 * t184;
t209 = qJD(2) * t135;
t207 = qJD(2) * t160;
t206 = qJD(4) * t152;
t205 = t132 * qJD(2);
t127 = t168 * qJD(5);
t148 = -t158 ^ 2 + t160 ^ 2;
t204 = t148 * qJD(2);
t203 = t158 * qJD(3);
t202 = t159 * qJD(2);
t201 = t160 * qJD(3);
t200 = t161 * qJD(2);
t154 = qJD(4) * t195;
t199 = t154 + qJD(5);
t194 = pkin(2) * t158 * qJD(2);
t193 = pkin(2) * t207;
t192 = t157 * t242;
t191 = t83 * t251;
t188 = t135 * t202;
t94 = t168 * t209;
t186 = t135 * t208;
t185 = t158 * t207;
t82 = t155 * t135;
t163 = (-t164 / 0.2e1 + t256) * t117 + t73 * t251;
t169 = t151 * t181 + t254 * t270;
t2 = t163 + t169;
t7 = t70 * t73;
t177 = t2 * qJD(1) + t7 * qJD(2);
t172 = pkin(4) * t255 + t257 * t270;
t4 = t191 + t172;
t8 = t70 * t83;
t176 = t4 * qJD(1) + t8 * qJD(2);
t173 = pkin(4) * t256 - t180 * t257;
t61 = t181 + t210;
t170 = t61 * qJD(1) + t70 * t209;
t109 = (t250 * t149 + t151 * t157) * pkin(3);
t12 = t171 * t117 + t167 * t119;
t162 = (t164 * t196 + t180 * t252) * pkin(3) + t180 * t254 + t164 * t253;
t6 = t162 + t173;
t165 = t12 * qJD(1) + t6 * qJD(2) + t109 * qJD(3);
t156 = qJ(5) * qJD(5);
t153 = qJD(3) * t195;
t147 = t155 * qJ(5);
t146 = t149 * qJD(5);
t81 = t155 * t168;
t69 = qJD(3) * t247;
t68 = t255 + t210;
t67 = t166 * t161 + t181;
t66 = t182 + t212;
t65 = t183 + t211;
t54 = t64 * qJD(2);
t52 = t63 * qJD(2);
t50 = t62 * qJD(2);
t45 = t155 * t247;
t31 = t155 * t62;
t30 = t67 * qJD(2) - t265;
t29 = t66 * qJD(2) + t266;
t28 = t65 * qJD(2) - t266;
t24 = t155 * t67 - t168 * t202;
t20 = -t230 / 0.2e1 + t229 / 0.2e1 - t231 / 0.2e1 - t248 / 0.2e1 + (t135 * t252 + t168 * t196) * pkin(3);
t11 = -t233 / 0.2e1 + t232 / 0.2e1 - t234 / 0.2e1 - t249 / 0.2e1 + (t117 * t252 + t119 * t196) * pkin(3);
t5 = t162 - t173;
t3 = t191 - t172;
t1 = t163 - t169;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t27; 0, 0, -t202, -t200, 0, 0, 0, 0, 0, -t160 * t202 - t161 * t203, t158 * t202 - t161 * t201, 0, 0, 0, 0, 0, t24, t155 * t66 + t188, t24, (t118 * t135 + t168 * t270) * qJD(2), t155 * t65 - t188, t219 + (t118 * t180 + t159 * t70 + t164 * t270) * qJD(2) + t1 * qJD(3) + t3 * qJD(4) + t68 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158 * t200 - t159 * t201, t159 * t203 - t160 * t200, 0, 0, 0, 0, 0, t30, t29, t30, 0, t28, t1 * qJD(2) + (t232 - t233) * qJD(3) + t11 * qJD(4) + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t29, t30, 0, t28, t3 * qJD(2) + t11 * qJD(3) + (-t234 - t249) * qJD(4) + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * qJD(2) + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t155 * t64, -t31, 0, -t155 * t63, qJD(3) * t2 + qJD(4) * t4 - qJD(5) * t61 - t219; 0, 0, 0, 0, t158 * t201, t148 * qJD(3), 0, 0, 0, -pkin(2) * t203, -pkin(2) * t201, t168 * t82, t155 * t46, 0, 0, 0, t59 * qJD(3) + t135 * t206, t60 * qJD(3) + t168 * t206, t22 * qJD(3) + t25 * qJD(4) + t135 * t127, 0, t23 * qJD(3) + t26 * qJD(4) + t132 * qJD(5), t7 * qJD(3) + t8 * qJD(4) - qJD(5) * t240; 0, 0, 0, 0, t185, t204, t201, -t203, 0, -pkin(6) * t201 - t194, pkin(6) * t203 - t193, t94, t218, t81, -t82, 0, t269 + t217, t268 - t263, t269 + t223, (t229 - t230) * qJD(3) + t20 * qJD(4) + t127, -t268 - t261, (-t149 * t180 + t151 * t164) * qJD(3) + t5 * qJD(4) + t259 + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t218, t81, -t82, 0, t186 + t269, t268 - t264, t269 + t221, t20 * qJD(3) + qJD(4) * t178 + t127, -t268 - t262, t5 * qJD(3) + (-pkin(4) * t164 - qJ(5) * t180) * qJD(4) + t259 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t81, t205, -t170 + t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t54, t50, 0, t52, -qJD(2) * t2 + qJD(4) * t12; 0, 0, 0, 0, -t185, -t204, 0, 0, 0, t194, t193, -t94, -t218, 0, 0, 0, -t217 + t215, t263, -t223 + t215, qJD(4) * t21, t261, t6 * qJD(4) - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, -t154, -t192, 0, t199, t109 * qJD(4) + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t153 - t154, -t45, t224, t199 + t153, (-pkin(4) * t157 + t250 * qJ(5)) * t242 + t146 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t149 * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t54, t50, 0, t52, -qJD(2) * t4 - qJD(3) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t218, 0, 0, 0, -t186 + t215, t264, -t221 + t215, -qJD(3) * t21, t262, -qJD(3) * t6 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t153, t69, -t224, qJD(5) - t153, t156 - t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, 0, -t205, t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -qJ(5) * qJD(4) - t149 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
