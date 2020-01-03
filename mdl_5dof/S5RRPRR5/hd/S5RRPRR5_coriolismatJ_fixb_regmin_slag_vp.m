% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:04:00
% EndTime: 2020-01-03 12:04:06
% DurationCPUTime: 2.10s
% Computational Cost: add. (2679->188), mult. (5493->239), div. (0->0), fcn. (6199->8), ass. (0->172)
t302 = qJD(4) + qJD(5);
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t164 = sin(pkin(9));
t263 = cos(qJ(4));
t210 = t263 * t164;
t165 = cos(pkin(9));
t167 = sin(qJ(4));
t237 = t167 * t165;
t141 = t210 + t237;
t137 = t141 * pkin(8);
t150 = (-pkin(7) - qJ(3)) * t164;
t161 = t165 * pkin(7);
t151 = t165 * qJ(3) + t161;
t196 = -t263 * t150 + t167 * t151;
t271 = -t137 - t196;
t179 = -t167 * t150 - t263 * t151;
t209 = t263 * t165;
t238 = t167 * t164;
t178 = t209 - t238;
t260 = t178 * pkin(8);
t81 = -t179 + t260;
t306 = t302 * (-t166 * t271 - t169 * t81);
t305 = t302 * (t166 * t81 - t169 * t271);
t168 = sin(qJ(2));
t259 = t168 * pkin(1);
t158 = qJ(3) + t259;
t135 = (-pkin(7) - t158) * t164;
t136 = t165 * t158 + t161;
t197 = -t263 * t135 + t167 * t136;
t270 = -t137 - t197;
t180 = -t167 * t135 - t263 * t136;
t67 = -t180 + t260;
t304 = t302 * (-t166 * t270 - t169 * t67);
t303 = t302 * (t166 * t67 - t169 * t270);
t217 = qJD(1) + qJD(2);
t100 = t166 * t141 - t169 * t178;
t128 = t169 * t141;
t239 = t166 * t178;
t272 = t128 + t239;
t284 = t100 ^ 2 - t272 ^ 2;
t293 = t217 * t284;
t292 = -t169 / 0.2e1;
t287 = t217 * t100;
t277 = t217 * t272;
t223 = t100 * qJD(5);
t286 = -t100 * qJD(4) - t223;
t285 = t100 * qJD(3);
t204 = t128 / 0.2e1;
t69 = 0.2e1 * t204 + t239;
t280 = t217 * t69;
t82 = -t141 ^ 2 + t178 ^ 2;
t279 = t217 * t82;
t98 = t204 - t128 / 0.2e1;
t278 = t217 * t98;
t159 = -t165 * pkin(3) - pkin(2);
t170 = cos(qJ(2));
t258 = t170 * pkin(1);
t149 = t159 - t258;
t202 = t159 / 0.2e1 + t149 / 0.2e1;
t276 = t202 * t141;
t114 = -pkin(4) * t178 + t159;
t107 = t114 - t258;
t203 = t114 / 0.2e1 + t107 / 0.2e1;
t275 = t203 * t272;
t274 = t217 * t178;
t195 = t217 * t141;
t162 = t164 ^ 2;
t163 = t165 ^ 2;
t154 = t162 + t163;
t273 = t217 * t154;
t117 = t141 * t258;
t118 = t178 * t258;
t229 = t118 * t292 + t166 * t117 / 0.2e1;
t230 = -t166 * t118 / 0.2e1 + t117 * t292;
t262 = pkin(2) * t168;
t261 = pkin(4) * t141;
t257 = pkin(1) * qJD(1);
t256 = pkin(1) * qJD(2);
t255 = pkin(4) * qJD(5);
t254 = qJD(4) * pkin(4);
t216 = t100 * t261;
t46 = t107 * t272;
t26 = -t46 - t216;
t234 = t26 * qJD(1);
t215 = t272 * t261;
t47 = t107 * t100;
t27 = t47 - t215;
t233 = t27 * qJD(1);
t119 = t154 * t158;
t94 = (-t262 + (t119 - t259) * t170) * pkin(1);
t232 = t94 * qJD(1);
t231 = t98 * qJD(5);
t198 = t154 * t170;
t134 = pkin(1) * t198;
t152 = t154 * qJD(3);
t228 = t134 * qJD(2) + t152;
t227 = qJD(1) * t107;
t226 = qJD(1) * t149;
t225 = qJD(2) * t114;
t224 = qJD(2) * t159;
t221 = t272 * qJD(5);
t220 = t119 * qJD(1);
t219 = t134 * qJD(1);
t130 = t178 * qJD(4);
t218 = t141 * qJD(4);
t214 = t168 * t256;
t213 = t168 * t257;
t208 = t100 * t227;
t207 = t272 * t227;
t206 = t178 * t226;
t205 = t141 * t226;
t200 = pkin(1) * t217;
t199 = pkin(4) * t302;
t193 = t100 * t213;
t192 = t272 * t213;
t191 = t178 * t213;
t190 = t141 * t213;
t189 = t164 * t213;
t188 = t168 * t200;
t61 = t114 * t272;
t28 = -t61 - t216;
t182 = t216 + t46 / 0.2e1 + t61 / 0.2e1;
t8 = t182 - t230;
t185 = t8 * qJD(1) - t28 * qJD(2);
t62 = t114 * t100;
t181 = t215 - t47 / 0.2e1 - t62 / 0.2e1;
t10 = t181 - t229;
t29 = t62 - t215;
t184 = t10 * qJD(1) - t29 * qJD(2);
t33 = qJD(4) * t272 + t69 * qJD(5);
t148 = t154 * qJ(3);
t171 = (qJ(3) + t158) * (t162 / 0.2e1 + t163 / 0.2e1);
t86 = -t259 / 0.2e1 + t171;
t183 = t86 * qJD(1) + t148 * qJD(2);
t19 = t203 * t100 + t229;
t177 = -t19 * qJD(1) - t100 * t225;
t18 = t230 - t275;
t176 = -t18 * qJD(1) + t225 * t272;
t172 = (-t209 / 0.2e1 + t238 / 0.2e1) * t258;
t39 = -t178 * t202 + t172;
t175 = -t39 * qJD(1) + t178 * t224;
t173 = (-t237 / 0.2e1 - t210 / 0.2e1) * t258;
t38 = t173 - t276;
t174 = -t38 * qJD(1) + t141 * t224;
t153 = t164 * t214;
t133 = t141 * qJD(3);
t129 = t178 * qJD(3);
t121 = t141 * t214;
t120 = t178 * t214;
t104 = t178 * t218;
t95 = t272 * qJD(3);
t92 = t98 * qJD(3);
t91 = t98 * qJD(4);
t85 = t259 / 0.2e1 + t171;
t84 = t272 * t214;
t83 = t100 * t214;
t78 = t82 * qJD(4);
t64 = t69 * qJD(3);
t48 = t178 * t195;
t41 = t173 + t276;
t40 = t172 + (t149 + t159) * t178 / 0.2e1;
t32 = -t69 * qJD(4) - t221;
t24 = t272 * t287;
t23 = t286 * t272;
t22 = t100 * t277;
t21 = t230 + t275;
t20 = t229 - (t107 + t114) * t100 / 0.2e1;
t11 = t181 + t229;
t9 = t182 + t230;
t3 = t302 * t284;
t1 = [0, 0, 0, 0, -t214, -t170 * t256, -t165 * t214, t153, t228, t94 * qJD(2) + t119 * qJD(3), t104, t78, 0, 0, 0, t149 * t218 - t120, t130 * t149 + t121, t23, t3, 0, 0, 0, -t26 * qJD(4) + t107 * t221 + t83, -t27 * qJD(4) - t107 * t223 + t84; 0, 0, 0, 0, -t188, -t170 * t200, -t165 * t188, t153 + t189, t219 + t228, t232 + t85 * qJD(3) + (qJ(3) * t198 - t262) * t256, t104, t78, 0, 0, 0, t41 * qJD(4) - t120 - t191, t40 * qJD(4) + t121 + t190, t23, t3, 0, 0, 0, t9 * qJD(4) + t21 * qJD(5) + t193 + t83, t11 * qJD(4) + t20 * qJD(5) + t192 + t84; 0, 0, 0, 0, 0, 0, 0, 0, t273, t85 * qJD(2) + t220, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t279, t130, -t218, 0, t41 * qJD(2) + qJD(4) * t180 + t205, t40 * qJD(2) + qJD(4) * t197 + t206, -t24, t293, t286, -t33, 0, t9 * qJD(2) - t234 + t304, t11 * qJD(2) - t233 + t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t293, t286, t32, 0, t21 * qJD(2) + t207 + t304 + t92, t20 * qJD(2) - t208 + t303; 0, 0, 0, 0, t213, t170 * t257, t165 * t213, -t189, t152 - t219, t86 * qJD(3) - t232, t104, t78, 0, 0, 0, -t38 * qJD(4) + t191, -t39 * qJD(4) - t190, t23, t3, 0, 0, 0, t8 * qJD(4) - t18 * qJD(5) - t193, t10 * qJD(4) - t19 * qJD(5) - t192; 0, 0, 0, 0, 0, 0, 0, 0, t152, t148 * qJD(3), t104, t78, 0, 0, 0, t159 * t218, t159 * t130, t23, t3, 0, 0, 0, -t28 * qJD(4) + t114 * t221, -t29 * qJD(4) - t114 * t223; 0, 0, 0, 0, 0, 0, 0, 0, t273, t183, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t279, t130, -t218, 0, qJD(4) * t179 + t174, qJD(4) * t196 + t175, -t24, t293, t286, -t33, 0, t185 + t306, t184 + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t293, t286, t32, 0, t176 + t92 + t306, t177 + t305; 0, 0, 0, 0, 0, 0, 0, 0, -t273, -t86 * qJD(2) - t220, 0, 0, 0, 0, 0, t218, t130, 0, 0, 0, 0, 0, t33, t286; 0, 0, 0, 0, 0, 0, 0, 0, -t273, -t183, 0, 0, 0, 0, 0, t218, t130, 0, 0, 0, 0, 0, t33, t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t274, 0, 0, 0, 0, 0, t277, -t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, -t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t279, 0, 0, 0, t38 * qJD(2) - t133 - t205, t39 * qJD(2) - t129 - t206, t24, -t293, 0, -t231, 0, -t8 * qJD(2) + t234 - t95, -t10 * qJD(2) + t233 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t279, 0, 0, 0, -t133 - t174, -t129 - t175, t24, -t293, 0, -t231, 0, -t185 - t95, -t184 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, -t274, 0, 0, 0, 0, 0, -t277, t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166 * t255, -t169 * t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, 0, -t166 * t199, -t169 * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t293, 0, t91, 0, t18 * qJD(2) - t207 - t64, t19 * qJD(2) + t208 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t293, 0, t91, 0, -t176 - t64, -t177 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, 0, t166 * t254, t169 * t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
