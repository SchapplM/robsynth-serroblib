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
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:02:59
% EndTime: 2022-01-20 11:03:04
% DurationCPUTime: 2.29s
% Computational Cost: add. (2677->186), mult. (5481->237), div. (0->0), fcn. (6191->8), ass. (0->170)
t300 = qJD(4) + qJD(5);
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t163 = sin(pkin(9));
t261 = cos(qJ(4));
t208 = t261 * t163;
t164 = cos(pkin(9));
t166 = sin(qJ(4));
t235 = t166 * t164;
t141 = t208 + t235;
t137 = t141 * pkin(8);
t150 = (-pkin(7) - qJ(3)) * t163;
t160 = t164 * pkin(7);
t151 = t164 * qJ(3) + t160;
t194 = -t261 * t150 + t166 * t151;
t269 = -t137 - t194;
t178 = -t166 * t150 - t261 * t151;
t207 = t261 * t164;
t236 = t166 * t163;
t177 = t207 - t236;
t258 = t177 * pkin(8);
t81 = -t178 + t258;
t304 = t300 * (-t165 * t269 - t168 * t81);
t303 = t300 * (t165 * t81 - t168 * t269);
t167 = sin(qJ(2));
t257 = t167 * pkin(1);
t157 = qJ(3) + t257;
t135 = (-pkin(7) - t157) * t163;
t136 = t164 * t157 + t160;
t195 = -t261 * t135 + t166 * t136;
t268 = -t137 - t195;
t179 = -t166 * t135 - t261 * t136;
t67 = -t179 + t258;
t302 = t300 * (-t165 * t268 - t168 * t67);
t301 = t300 * (t165 * t67 - t168 * t268);
t215 = qJD(1) + qJD(2);
t100 = t165 * t141 - t168 * t177;
t128 = t168 * t141;
t237 = t165 * t177;
t270 = t128 + t237;
t282 = t100 ^ 2 - t270 ^ 2;
t291 = t215 * t282;
t290 = -t168 / 0.2e1;
t285 = t215 * t100;
t275 = t215 * t270;
t221 = t100 * qJD(5);
t284 = -t100 * qJD(4) - t221;
t283 = t100 * qJD(3);
t202 = t128 / 0.2e1;
t69 = 0.2e1 * t202 + t237;
t278 = t215 * t69;
t82 = -t141 ^ 2 + t177 ^ 2;
t277 = t215 * t82;
t98 = t202 - t128 / 0.2e1;
t276 = t215 * t98;
t158 = -t164 * pkin(3) - pkin(2);
t169 = cos(qJ(2));
t256 = t169 * pkin(1);
t149 = t158 - t256;
t200 = t158 / 0.2e1 + t149 / 0.2e1;
t274 = t200 * t141;
t114 = -pkin(4) * t177 + t158;
t107 = t114 - t256;
t201 = t114 / 0.2e1 + t107 / 0.2e1;
t273 = t201 * t270;
t272 = t215 * t177;
t193 = t215 * t141;
t161 = t163 ^ 2;
t162 = t164 ^ 2;
t153 = t161 + t162;
t271 = t215 * t153;
t117 = t141 * t256;
t118 = t177 * t256;
t227 = t118 * t290 + t165 * t117 / 0.2e1;
t228 = -t165 * t118 / 0.2e1 + t117 * t290;
t260 = pkin(2) * t167;
t259 = pkin(4) * t141;
t255 = pkin(1) * qJD(1);
t254 = pkin(1) * qJD(2);
t253 = pkin(4) * qJD(5);
t252 = qJD(4) * pkin(4);
t214 = t100 * t259;
t46 = t107 * t270;
t26 = -t46 - t214;
t232 = t26 * qJD(1);
t213 = t270 * t259;
t47 = t107 * t100;
t27 = t47 - t213;
t231 = t27 * qJD(1);
t119 = t153 * t157;
t94 = (-t260 + (t119 - t257) * t169) * pkin(1);
t230 = t94 * qJD(1);
t229 = t98 * qJD(5);
t196 = t153 * t169;
t134 = pkin(1) * t196;
t152 = t153 * qJD(3);
t226 = t134 * qJD(2) + t152;
t225 = qJD(1) * t107;
t224 = qJD(1) * t149;
t223 = qJD(2) * t114;
t222 = qJD(2) * t158;
t219 = t270 * qJD(5);
t218 = t119 * qJD(1);
t217 = t134 * qJD(1);
t130 = t177 * qJD(4);
t216 = t141 * qJD(4);
t212 = t167 * t254;
t211 = t167 * t255;
t206 = t100 * t225;
t205 = t270 * t225;
t204 = t177 * t224;
t203 = t141 * t224;
t198 = pkin(1) * t215;
t197 = pkin(4) * t300;
t191 = t100 * t211;
t190 = t270 * t211;
t189 = t177 * t211;
t188 = t141 * t211;
t187 = t167 * t198;
t61 = t114 * t270;
t28 = -t61 - t214;
t181 = t214 + t46 / 0.2e1 + t61 / 0.2e1;
t8 = t181 - t228;
t184 = t8 * qJD(1) - t28 * qJD(2);
t62 = t114 * t100;
t180 = t213 - t47 / 0.2e1 - t62 / 0.2e1;
t10 = t180 - t227;
t29 = t62 - t213;
t183 = t10 * qJD(1) - t29 * qJD(2);
t33 = qJD(4) * t270 + t69 * qJD(5);
t148 = t153 * qJ(3);
t170 = (qJ(3) + t157) * (t161 / 0.2e1 + t162 / 0.2e1);
t86 = -t257 / 0.2e1 + t170;
t182 = t86 * qJD(1) + t148 * qJD(2);
t19 = t201 * t100 + t227;
t176 = -t19 * qJD(1) - t100 * t223;
t18 = t228 - t273;
t175 = -t18 * qJD(1) + t223 * t270;
t171 = (-t207 / 0.2e1 + t236 / 0.2e1) * t256;
t39 = -t177 * t200 + t171;
t174 = -t39 * qJD(1) + t177 * t222;
t172 = (-t235 / 0.2e1 - t208 / 0.2e1) * t256;
t38 = t172 - t274;
t173 = -t38 * qJD(1) + t141 * t222;
t133 = t141 * qJD(3);
t129 = t177 * qJD(3);
t121 = t141 * t212;
t120 = t177 * t212;
t104 = t177 * t216;
t95 = t270 * qJD(3);
t92 = t98 * qJD(3);
t91 = t98 * qJD(4);
t85 = t257 / 0.2e1 + t170;
t84 = t270 * t212;
t83 = t100 * t212;
t78 = t82 * qJD(4);
t64 = t69 * qJD(3);
t48 = t177 * t193;
t41 = t172 + t274;
t40 = t171 + (t149 + t158) * t177 / 0.2e1;
t32 = -t69 * qJD(4) - t219;
t24 = t270 * t285;
t23 = t284 * t270;
t22 = t100 * t275;
t21 = t228 + t273;
t20 = t227 - (t107 + t114) * t100 / 0.2e1;
t11 = t180 + t227;
t9 = t181 + t228;
t3 = t300 * t282;
t1 = [0, 0, 0, 0, -t212, -t169 * t254, -t164 * t212, t226, t94 * qJD(2) + t119 * qJD(3), t104, t78, 0, 0, 0, t149 * t216 - t120, t149 * t130 + t121, t23, t3, 0, 0, 0, -t26 * qJD(4) + t107 * t219 + t83, -t27 * qJD(4) - t107 * t221 + t84; 0, 0, 0, 0, -t187, -t169 * t198, -t164 * t187, t217 + t226, t230 + t85 * qJD(3) + (qJ(3) * t196 - t260) * t254, t104, t78, 0, 0, 0, t41 * qJD(4) - t120 - t189, t40 * qJD(4) + t121 + t188, t23, t3, 0, 0, 0, t9 * qJD(4) + t21 * qJD(5) + t191 + t83, t11 * qJD(4) + t20 * qJD(5) + t190 + t84; 0, 0, 0, 0, 0, 0, 0, t271, t85 * qJD(2) + t218, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t277, t130, -t216, 0, t41 * qJD(2) + t179 * qJD(4) + t203, t40 * qJD(2) + t195 * qJD(4) + t204, -t24, t291, t284, -t33, 0, t9 * qJD(2) - t232 + t302, t11 * qJD(2) - t231 + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t291, t284, t32, 0, t21 * qJD(2) + t205 + t302 + t92, t20 * qJD(2) - t206 + t301; 0, 0, 0, 0, t211, t169 * t255, t164 * t211, t152 - t217, t86 * qJD(3) - t230, t104, t78, 0, 0, 0, -t38 * qJD(4) + t189, -t39 * qJD(4) - t188, t23, t3, 0, 0, 0, t8 * qJD(4) - t18 * qJD(5) - t191, t10 * qJD(4) - t19 * qJD(5) - t190; 0, 0, 0, 0, 0, 0, 0, t152, t148 * qJD(3), t104, t78, 0, 0, 0, t158 * t216, t158 * t130, t23, t3, 0, 0, 0, -t28 * qJD(4) + t114 * t219, -t29 * qJD(4) - t114 * t221; 0, 0, 0, 0, 0, 0, 0, t271, t182, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t277, t130, -t216, 0, t178 * qJD(4) + t173, t194 * qJD(4) + t174, -t24, t291, t284, -t33, 0, t184 + t304, t183 + t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t291, t284, t32, 0, t175 + t92 + t304, t176 + t303; 0, 0, 0, 0, 0, 0, 0, -t271, -t86 * qJD(2) - t218, 0, 0, 0, 0, 0, t216, t130, 0, 0, 0, 0, 0, t33, t284; 0, 0, 0, 0, 0, 0, 0, -t271, -t182, 0, 0, 0, 0, 0, t216, t130, 0, 0, 0, 0, 0, t33, t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t272, 0, 0, 0, 0, 0, t275, -t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, -t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t277, 0, 0, 0, t38 * qJD(2) - t133 - t203, t39 * qJD(2) - t129 - t204, t24, -t291, 0, -t229, 0, -t8 * qJD(2) + t232 - t95, -t10 * qJD(2) + t231 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t277, 0, 0, 0, -t133 - t173, -t129 - t174, t24, -t291, 0, -t229, 0, -t184 - t95, -t183 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, -t272, 0, 0, 0, 0, 0, -t275, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 * t253, -t168 * t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t276, 0, -t165 * t197, -t168 * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t291, 0, t91, 0, t18 * qJD(2) - t205 - t64, t19 * qJD(2) + t206 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t291, 0, t91, 0, -t175 - t64, -t176 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, 0, t165 * t252, t168 * t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
