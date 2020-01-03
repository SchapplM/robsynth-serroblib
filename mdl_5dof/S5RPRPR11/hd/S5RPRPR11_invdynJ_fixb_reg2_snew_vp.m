% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:28:04
% DurationCPUTime: 3.61s
% Computational Cost: add. (6965->325), mult. (17292->416), div. (0->0), fcn. (12338->8), ass. (0->189)
t172 = sin(pkin(8));
t173 = cos(pkin(8));
t178 = cos(qJ(3));
t175 = sin(qJ(3));
t208 = qJD(1) * t178;
t209 = qJD(1) * t175;
t148 = t172 * t209 - t173 * t208;
t150 = t172 * t208 + t173 * t209;
t215 = t150 * t148;
t242 = qJDD(3) + t215;
t222 = t242 * t175;
t145 = t150 ^ 2;
t179 = qJD(3) ^ 2;
t246 = -t145 - t179;
t74 = -t178 * t246 + t222;
t221 = t242 * t178;
t76 = t175 * t246 + t221;
t277 = qJ(2) * (t172 * t74 - t173 * t76);
t276 = pkin(6) * t74;
t275 = pkin(6) * t76;
t240 = t148 ^ 2;
t131 = t240 - t179;
t273 = t172 * (-t131 * t178 + t222) - t173 * (t131 * t175 + t221);
t243 = qJDD(3) - t215;
t220 = t243 * t175;
t244 = -t179 - t240;
t252 = t178 * t244 - t220;
t270 = pkin(6) * t252;
t203 = t173 * qJDD(1);
t204 = t172 * qJDD(1);
t146 = t175 * t204 - t178 * t203;
t211 = t172 * t178;
t147 = (t173 * t175 + t211) * qJDD(1);
t254 = -t146 * t175 - t178 * t147;
t269 = pkin(6) * t254;
t104 = t178 * t243;
t255 = t175 * t244 + t104;
t268 = pkin(6) * t255;
t180 = qJD(1) ^ 2;
t176 = sin(qJ(1));
t237 = cos(qJ(1));
t188 = t237 * g(1) + t176 * g(2);
t256 = -t180 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t188;
t253 = -t146 * t178 + t147 * t175;
t86 = -t240 - t145;
t265 = -pkin(2) * t86 + pkin(6) * t253;
t207 = t150 * qJD(3);
t118 = t146 + 0.2e1 * t207;
t264 = qJ(2) * (-t172 * t255 + t173 * t252) - pkin(1) * t118;
t132 = -t145 + t179;
t263 = t172 * (-t132 * t175 + t104) + t173 * (t178 * t132 + t220);
t262 = qJ(2) * (-t172 * t254 + t173 * t253) - pkin(1) * t86;
t174 = sin(qJ(5));
t164 = qJDD(3) - qJDD(5);
t177 = cos(qJ(5));
t101 = -t177 * t148 + t150 * t174;
t103 = t148 * t174 + t150 * t177;
t72 = t103 * t101;
t250 = -t72 - t164;
t261 = t174 * t250;
t259 = t177 * t250;
t141 = qJD(3) * t148;
t121 = t147 - t141;
t248 = -t141 + t121;
t257 = t248 * qJ(4);
t119 = t146 + t207;
t61 = -qJD(5) * t101 + t119 * t174 + t121 * t177;
t167 = qJD(3) - qJD(5);
t95 = t101 * t167;
t251 = t61 + t95;
t197 = g(1) * t176 - t237 * g(2);
t190 = -qJDD(2) + t197;
t165 = t172 ^ 2;
t166 = t173 ^ 2;
t210 = t165 + t166;
t232 = t173 * pkin(2);
t112 = (pkin(1) + t232) * qJDD(1) + (t210 * pkin(6) + qJ(2)) * t180 + t190;
t249 = -pkin(3) * t207 + t112;
t247 = t141 + t121;
t245 = t145 - t240;
t225 = qJDD(1) * pkin(1);
t226 = qJ(2) * t180;
t143 = t190 + t225 + t226;
t214 = t166 * t180;
t241 = qJ(2) * t214 + t165 * t226 - t143 - t225;
t195 = -t177 * t119 + t121 * t174;
t45 = (qJD(5) + t167) * t103 + t195;
t99 = t101 ^ 2;
t100 = t103 ^ 2;
t163 = t167 ^ 2;
t239 = 2 * qJD(4);
t238 = pkin(3) + pkin(4);
t236 = pkin(3) * t119;
t235 = pkin(3) * t175;
t234 = pkin(3) * t178;
t233 = g(3) * t173;
t187 = -t233 + (-pkin(6) * qJDD(1) + t180 * t232 - t256) * t172;
t192 = -g(3) * t172 + t256 * t173;
t96 = -pkin(2) * t214 + pkin(6) * t203 + t192;
t64 = t175 * t187 + t178 * t96;
t63 = t175 * t96 - t178 * t187;
t35 = t175 * t64 - t178 * t63;
t231 = t172 * t35;
t128 = -qJD(3) * pkin(4) - pkin(7) * t150;
t184 = t249 + t257;
t27 = -t236 - t119 * pkin(4) - t240 * pkin(7) + (t239 + t128) * t150 + t184;
t230 = t174 * t27;
t68 = -t72 + t164;
t229 = t174 * t68;
t228 = t177 * t27;
t227 = t177 * t68;
t224 = t112 * t175;
t223 = t112 * t178;
t219 = t118 * t175;
t218 = t118 * t178;
t213 = t167 * t174;
t212 = t167 * t177;
t202 = t150 * t239;
t198 = -qJ(4) * t175 - pkin(2);
t107 = pkin(3) * t148 - qJ(4) * t150;
t51 = -qJDD(3) * pkin(3) - t179 * qJ(4) + t107 * t150 + qJDD(4) + t63;
t29 = -t243 * pkin(4) - t247 * pkin(7) + t51;
t191 = qJDD(3) * qJ(4) + qJD(3) * t239 - t148 * t107 + t64;
t50 = -pkin(3) * t179 + t191;
t30 = -pkin(4) * t240 + pkin(7) * t119 + qJD(3) * t128 + t50;
t15 = t174 * t30 - t177 * t29;
t36 = t175 * t63 + t178 * t64;
t194 = t172 * (t256 * t172 + t233) + t173 * t192;
t16 = t174 * t29 + t177 * t30;
t7 = -t15 * t177 + t16 * t174;
t8 = t15 * t174 + t16 * t177;
t186 = t172 * (t119 * t175 + t178 * t141) + t173 * (-t178 * t119 + t175 * t141);
t129 = t175 * t207;
t185 = t172 * t129 + (-t148 * t211 + t173 * (-t148 * t175 - t150 * t178)) * qJD(3);
t183 = t184 + t202;
t160 = t166 * qJDD(1);
t159 = t165 * qJDD(1);
t152 = t210 * t180;
t120 = t147 - 0.2e1 * t141;
t93 = -t100 + t163;
t92 = t99 - t163;
t90 = -t100 - t163;
t71 = t100 - t99;
t67 = -t163 - t99;
t66 = (t101 * t177 - t103 * t174) * t167;
t65 = (-t101 * t174 - t103 * t177) * t167;
t60 = -qJD(5) * t103 - t195;
t59 = -t99 - t100;
t58 = t172 * (t121 * t178 - t129) + t173 * (t121 * t175 + t178 * t207);
t57 = t177 * t92 + t229;
t56 = -t174 * t93 + t259;
t55 = -t174 * t92 + t227;
t54 = -t177 * t93 - t261;
t53 = -t174 * t90 + t227;
t52 = t177 * t90 + t229;
t49 = t61 - t95;
t44 = (qJD(5) - t167) * t103 + t195;
t43 = t183 - t236;
t42 = t103 * t213 + t177 * t61;
t41 = t103 * t212 - t174 * t61;
t40 = -t101 * t212 - t174 * t60;
t39 = t101 * t213 - t177 * t60;
t38 = t177 * t67 - t261;
t37 = t174 * t67 + t259;
t34 = (-t118 - t119) * pkin(3) + t183;
t33 = -qJ(4) * t86 + t51;
t32 = (-t179 - t86) * pkin(3) + t191;
t31 = t202 - t236 + t249 + 0.2e1 * t257;
t26 = t175 * t52 + t178 * t53;
t25 = t175 * t53 - t178 * t52;
t22 = t174 * t49 - t177 * t45;
t21 = -t174 * t251 - t177 * t44;
t20 = -t174 * t45 - t177 * t49;
t19 = t174 * t44 - t177 * t251;
t18 = t175 * t37 + t178 * t38;
t17 = t175 * t38 - t178 * t37;
t14 = -pkin(7) * t52 + qJ(4) * t251 + t228;
t13 = -pkin(7) * t37 + qJ(4) * t44 + t230;
t12 = t175 * t20 + t178 * t22;
t11 = t175 * t22 - t178 * t20;
t10 = -pkin(7) * t53 + t238 * t251 - t230;
t9 = -pkin(7) * t38 + t238 * t44 + t228;
t6 = -pkin(7) * t7 + qJ(4) * t27;
t5 = -pkin(7) * t20 + qJ(4) * t59 - t7;
t4 = -pkin(7) * t8 + t238 * t27;
t3 = -pkin(7) * t22 + t238 * t59 - t8;
t2 = t175 * t7 + t178 * t8;
t1 = t175 * t8 - t178 * t7;
t23 = [0, 0, 0, 0, 0, qJDD(1), t197, t188, 0, 0, t159, 0.2e1 * t172 * t203, 0, t160, 0, 0, -t241 * t173, t241 * t172, pkin(1) * t152 + qJ(2) * (t160 + t159) + t194, pkin(1) * t143 + qJ(2) * t194, t58, t172 * (-t120 * t175 - t218) + t173 * (t120 * t178 - t219), t263, t186, -t273, t185, t172 * (-t224 - t268) + t173 * (-pkin(2) * t118 + t223 + t270) + t264, t172 * (-t223 + t276) + t173 * (-pkin(2) * t120 - t224 - t275) - pkin(1) * t120 + t277, t172 * (-t35 - t269) + t173 * (t265 + t36) + t262, -pkin(6) * t231 + t173 * (pkin(2) * t112 + pkin(6) * t36) + pkin(1) * t112 + qJ(2) * (t173 * t36 - t231), t58, t263, t172 * (t175 * t248 + t218) + t173 * (-t178 * t248 + t219), t185, t273, t186, t172 * (-qJ(4) * t218 - t175 * t34 - t268) + t173 * (t198 * t118 + t178 * t34 + t270) + t264, t172 * (-t175 * t32 + t178 * t33 - t269) + t173 * (t175 * t33 + t178 * t32 + t265) + t262, t172 * (t178 * t31 - t276) + t173 * (t175 * t31 + t275) - t277 + (-t172 * t235 + t173 * (pkin(2) + t234) + pkin(1)) * t248, (t172 * (qJ(4) * t178 - t235) + t173 * (-t198 + t234) + pkin(1)) * t43 + (qJ(2) + pkin(6)) * (-t172 * (t175 * t50 - t178 * t51) + t173 * (t175 * t51 + t178 * t50)), t172 * (-t175 * t41 + t178 * t42) + t173 * (t175 * t42 + t178 * t41), t172 * (-t175 * t19 + t178 * t21) + t173 * (t175 * t21 + t178 * t19), t172 * (-t175 * t54 + t178 * t56) + t173 * (t175 * t56 + t178 * t54), t172 * (-t175 * t39 + t178 * t40) + t173 * (t175 * t40 + t178 * t39), t172 * (-t175 * t55 + t178 * t57) + t173 * (t175 * t57 + t178 * t55), t172 * (-t175 * t65 + t178 * t66) + t173 * (t175 * t66 + t178 * t65), t172 * (-pkin(6) * t17 + t13 * t178 - t175 * t9) + t173 * (pkin(2) * t44 + pkin(6) * t18 + t13 * t175 + t178 * t9) + pkin(1) * t44 + qJ(2) * (-t17 * t172 + t173 * t18), t172 * (-pkin(6) * t25 - t10 * t175 + t14 * t178) + t173 * (pkin(2) * t251 + pkin(6) * t26 + t10 * t178 + t14 * t175) + pkin(1) * t251 + qJ(2) * (-t172 * t25 + t173 * t26), t172 * (-pkin(6) * t11 - t175 * t3 + t178 * t5) + t173 * (pkin(2) * t59 + pkin(6) * t12 + t175 * t5 + t178 * t3) + pkin(1) * t59 + qJ(2) * (-t11 * t172 + t12 * t173), t172 * (-pkin(6) * t1 - t175 * t4 + t178 * t6) + t173 * (pkin(2) * t27 + pkin(6) * t2 + t175 * t6 + t178 * t4) + pkin(1) * t27 + qJ(2) * (-t1 * t172 + t173 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, t204, -t152, -t143, 0, 0, 0, 0, 0, 0, t118, t120, t86, -t112, 0, 0, 0, 0, 0, 0, t118, t86, -t248, -t43, 0, 0, 0, 0, 0, 0, -t44, -t251, -t59, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, t245, t147, -t215, -t146, qJDD(3), -t63, -t64, 0, 0, t215, t247, -t245, qJDD(3), t146, -t215, pkin(3) * t243 + qJ(4) * t244 - t51, -pkin(3) * t147 - qJ(4) * t146, qJ(4) * t242 + (-t246 - t179) * pkin(3) + t191, -pkin(3) * t51 + qJ(4) * t50, -t72, -t71, -t49, t72, t45, t164, qJ(4) * t38 - t238 * t37 + t15, qJ(4) * t53 - t238 * t52 + t16, qJ(4) * t22 - t238 * t20, qJ(4) * t8 - t238 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t147, t246, t51, 0, 0, 0, 0, 0, 0, t37, t52, t20, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t71, t49, -t72, -t45, -t164, -t15, -t16, 0, 0;];
tauJ_reg = t23;
