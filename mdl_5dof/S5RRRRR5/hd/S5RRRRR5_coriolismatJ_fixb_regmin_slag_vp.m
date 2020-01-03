% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:25
% EndTime: 2020-01-03 12:13:32
% DurationCPUTime: 2.65s
% Computational Cost: add. (2636->311), mult. (5555->383), div. (0->0), fcn. (5200->8), ass. (0->243)
t242 = qJD(4) + qJD(5);
t187 = cos(qJ(3));
t188 = cos(qJ(2));
t300 = t188 * pkin(1);
t238 = pkin(2) + t300;
t161 = t187 * t238;
t184 = sin(qJ(3));
t185 = sin(qJ(2));
t280 = t184 * t185;
t134 = pkin(1) * t280 - t161;
t186 = cos(qJ(4));
t305 = cos(qJ(5));
t236 = t305 * t186;
t182 = sin(qJ(5));
t183 = sin(qJ(4));
t283 = t182 * t183;
t199 = t236 / 0.2e1 - t283 / 0.2e1;
t297 = t199 * t134;
t243 = -qJD(2) - qJD(3);
t220 = qJD(1) - t243;
t146 = -t236 + t283;
t237 = t305 * t183;
t282 = t182 * t186;
t148 = -t237 - t282;
t56 = t146 ^ 2 - t148 ^ 2;
t324 = t220 * t56;
t169 = -t183 ^ 2 + t186 ^ 2;
t323 = t220 * t169;
t315 = pkin(8) + pkin(9);
t159 = t315 * t183;
t160 = t315 * t186;
t322 = t242 * (t182 * t159 - t305 * t160);
t321 = t242 * (t305 * t159 + t182 * t160);
t95 = t184 * pkin(2);
t174 = pkin(8) + t95;
t298 = pkin(9) + t174;
t140 = t298 * t183;
t141 = t298 * t186;
t320 = t242 * (t182 * t140 - t305 * t141);
t319 = t242 * (t305 * t140 + t182 * t141);
t214 = t184 * t238;
t278 = t187 * t185;
t135 = pkin(1) * t278 + t214;
t131 = pkin(8) + t135;
t299 = pkin(9) + t131;
t86 = t299 * t183;
t87 = t299 * t186;
t318 = t242 * (t182 * t86 - t305 * t87);
t317 = t242 * (t182 * t87 + t305 * t86);
t316 = pkin(3) / 0.2e1;
t130 = -pkin(3) + t134;
t302 = t186 * pkin(4);
t103 = t130 - t302;
t314 = t103 / 0.2e1;
t313 = -t130 / 0.2e1;
t279 = t184 * t188;
t138 = (t278 + t279) * pkin(1);
t312 = -t138 / 0.2e1;
t176 = -pkin(3) - t302;
t301 = t187 * pkin(2);
t157 = t176 - t301;
t311 = t157 / 0.2e1;
t175 = -pkin(3) - t301;
t310 = -t175 / 0.2e1;
t309 = t176 / 0.2e1;
t308 = t183 / 0.2e1;
t307 = -t186 / 0.2e1;
t306 = t186 / 0.2e1;
t304 = pkin(4) * t182;
t303 = pkin(4) * t183;
t296 = pkin(1) * qJD(1);
t295 = pkin(1) * qJD(2);
t294 = pkin(2) * qJD(2);
t293 = pkin(2) * qJD(3);
t292 = pkin(3) * qJD(3);
t133 = t148 * t303;
t290 = t103 * t146;
t63 = -t290 / 0.2e1;
t291 = t63 - t133;
t289 = t103 * t148;
t288 = t138 * t186;
t287 = t157 * t146;
t286 = t157 * t148;
t285 = t176 * t146;
t284 = t176 * t148;
t277 = t187 * t188;
t139 = (t277 - t280) * pkin(1);
t281 = t183 * t139;
t132 = t146 * t303;
t40 = t132 - t289;
t268 = t40 * qJD(1);
t41 = -t133 - t290;
t267 = t41 * qJD(1);
t216 = t95 / 0.2e1 + t135 / 0.2e1;
t208 = t312 + t216;
t52 = t208 * t186;
t266 = t52 * qJD(1);
t261 = t95 * qJD(1);
t97 = t161 / 0.2e1 + (-t300 / 0.2e1 + pkin(2) / 0.2e1) * t187;
t260 = t97 * qJD(1);
t224 = t139 * t307;
t259 = t182 * t281 / 0.2e1 + t305 * t224;
t125 = t135 * qJD(3);
t136 = t138 * qJD(2);
t258 = -t136 - t125;
t257 = qJD(1) * t103;
t256 = qJD(1) * t130;
t255 = qJD(2) * t157;
t254 = qJD(2) * t175;
t253 = qJD(3) * t176;
t252 = qJD(5) * t103;
t251 = qJD(5) * t157;
t250 = qJD(5) * t176;
t249 = t134 * qJD(1);
t248 = t135 * qJD(1);
t247 = t138 * qJD(1);
t246 = t139 * qJD(1);
t245 = t183 * qJD(4);
t181 = t186 * qJD(4);
t244 = -qJD(1) - qJD(2);
t241 = t184 * t293;
t240 = t184 * t294;
t239 = -t301 / 0.2e1;
t235 = t146 * t257;
t234 = t148 * t257;
t233 = t183 * t256;
t232 = t186 * t256;
t231 = t146 * t248;
t230 = t148 * t248;
t229 = t183 * t248;
t228 = t146 * t247;
t227 = t148 * t247;
t226 = t183 * t247;
t223 = t309 + t311;
t222 = t305 * qJD(4);
t221 = t305 * qJD(5);
t219 = pkin(1) * t244;
t218 = pkin(2) * t243;
t113 = -t285 / 0.2e1;
t13 = t113 + t63 + t297;
t88 = -t287 / 0.2e1;
t17 = t88 + t63 + t259;
t75 = t242 * t148;
t215 = t134 / 0.2e1 + t316 + t313;
t213 = t184 * t218;
t212 = -t139 / 0.2e1 + t310 + t313;
t54 = t132 - t286;
t200 = t282 / 0.2e1 + t237 / 0.2e1;
t192 = t200 * t139;
t15 = (t311 + t314) * t148 - t192;
t7 = -t132 + t15;
t211 = -t7 * qJD(1) + t54 * qJD(2);
t55 = -t133 - t287;
t9 = t199 * t139 + t291 + t88;
t210 = t9 * qJD(1) + t55 * qJD(2);
t209 = t239 + t316 + t310;
t64 = t290 / 0.2e1;
t16 = t287 / 0.2e1 + t64 + t259;
t207 = -t16 * qJD(1) - t146 * t255;
t206 = -t15 * qJD(1) - t148 * t255;
t42 = t212 * t183;
t205 = t42 * qJD(1) - t183 * t254;
t43 = t212 * t186;
t204 = t43 * qJD(1) - t186 * t254;
t28 = t208 * t146;
t203 = -t28 * qJD(1) - t146 * t240;
t30 = t208 * t148;
t202 = t30 * qJD(1) + t148 * t240;
t51 = t208 * t183;
t201 = -t51 * qJD(1) - t183 * t240;
t190 = t200 * t301;
t31 = t223 * t148 - t190;
t23 = -t132 + t31;
t191 = t200 * t134;
t11 = (t309 + t314) * t148 + t191;
t3 = -t132 + t11;
t61 = t132 - t284;
t198 = t3 * qJD(1) + t23 * qJD(2) - t61 * qJD(3);
t189 = t199 * t301;
t32 = t223 * t146 - t189;
t24 = t133 + t32;
t5 = t113 + t291 - t297;
t62 = -t133 - t285;
t197 = t5 * qJD(1) - t24 * qJD(2) + t62 * qJD(3);
t46 = t215 * t183;
t94 = t209 * t183;
t196 = t46 * qJD(1) + t94 * qJD(2) + t183 * t292;
t47 = t215 * t186;
t96 = t209 * t186;
t195 = t47 * qJD(1) + t96 * qJD(2) + t186 * t292;
t194 = t11 * qJD(1) + t31 * qJD(2) + t148 * t253;
t12 = t285 / 0.2e1 + t64 + t297;
t193 = t12 * qJD(1) + t32 * qJD(2) + t146 * t253;
t65 = -t289 / 0.2e1;
t89 = -t286 / 0.2e1;
t18 = t89 + t65 - t192;
t114 = -t284 / 0.2e1;
t14 = t114 + t65 + t191;
t33 = t113 + t88 - t189;
t34 = t114 + t89 - t190;
t180 = pkin(3) * t307;
t179 = -pkin(3) * t183 / 0.2e1;
t170 = t183 * t181;
t168 = t183 * t241;
t158 = t169 * qJD(4);
t156 = t175 * t306;
t155 = t175 * t308;
t137 = t139 * qJD(2);
t129 = t148 * t241;
t128 = t146 * t241;
t124 = t134 * qJD(3);
t119 = t183 * t136;
t112 = t183 * t125;
t102 = t130 * t306;
t101 = t130 * t308;
t100 = t220 * t186 * t183;
t99 = t186 * t239 + t156 + t180;
t98 = t183 * t239 + t155 + t179;
t83 = t148 * t136;
t82 = t146 * t136;
t77 = t239 - t161 / 0.2e1 + (t280 - t277 / 0.2e1) * pkin(1);
t76 = -t95 / 0.2e1 - t214 / 0.2e1 + (-t278 - t279 / 0.2e1) * pkin(1);
t74 = t242 * t146;
t69 = t148 * t125;
t68 = t146 * t125;
t53 = -t288 / 0.2e1 - t216 * t186;
t50 = t138 * t308 + t216 * t183;
t49 = t134 * t306 + t102 + t180;
t48 = t134 * t308 + t101 + t179;
t45 = t156 + t102 + t224;
t44 = t155 + t101 - t281 / 0.2e1;
t37 = t146 * t75;
t29 = (-t216 + t312) * t148;
t27 = (t138 / 0.2e1 + t216) * t146;
t26 = -t133 + t33;
t25 = t132 + t34;
t22 = t220 * t148 * t146;
t21 = t242 * t56;
t10 = -t133 + t17;
t8 = t132 + t18;
t6 = -t133 + t13;
t4 = t132 + t14;
t1 = [0, 0, 0, 0, -t185 * t295, -t188 * t295, 0, t258, -t137 + t124, t170, t158, 0, 0, 0, t130 * t245 + t258 * t186, t130 * t181 + t112 + t119, t37, t21, 0, 0, 0, t40 * qJD(4) - t148 * t252 + t68 + t82, t41 * qJD(4) - t146 * t252 - t69 - t83; 0, 0, 0, 0, t185 * t219, t188 * t219, 0, t76 * qJD(3) - t136 - t247, t77 * qJD(3) - t137 - t246, t170, t158, 0, 0, 0, t53 * qJD(3) + t44 * qJD(4) + t244 * t288, t50 * qJD(3) + t45 * qJD(4) + t119 + t226, t37, t21, 0, 0, 0, t27 * qJD(3) + t8 * qJD(4) + t18 * qJD(5) + t228 + t82, t29 * qJD(3) + t10 * qJD(4) + t17 * qJD(5) - t227 - t83; 0, 0, 0, 0, 0, 0, 0, t76 * qJD(2) - t125 - t248, t77 * qJD(2) + t124 + t249, t170, t158, 0, 0, 0, t53 * qJD(2) + t48 * qJD(4) + (-qJD(1) - qJD(3)) * t186 * t135, t50 * qJD(2) + t49 * qJD(4) + t112 + t229, t37, t21, 0, 0, 0, t27 * qJD(2) + t4 * qJD(4) + t14 * qJD(5) + t231 + t68, t29 * qJD(2) + t6 * qJD(4) + t13 * qJD(5) - t230 - t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t323, t181, -t245, 0, t44 * qJD(2) + t48 * qJD(3) - t131 * t181 + t233, t45 * qJD(2) + t49 * qJD(3) + t131 * t245 + t232, t22, t324, -t74, t75, 0, t8 * qJD(2) + t4 * qJD(3) + t268 + t318, t10 * qJD(2) + t6 * qJD(3) + t267 + t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t324, -t74, t75, 0, t18 * qJD(2) + t14 * qJD(3) - t234 + t318, t17 * qJD(2) + t13 * qJD(3) - t235 + t317; 0, 0, 0, 0, t185 * t296, t188 * t296, 0, -t95 * qJD(3) + t247, -t97 * qJD(3) + t246, t170, t158, 0, 0, 0, -t52 * qJD(3) - t42 * qJD(4) + t186 * t247, t51 * qJD(3) - t43 * qJD(4) - t226, t37, t21, 0, 0, 0, t28 * qJD(3) - t7 * qJD(4) - t15 * qJD(5) - t228, -t30 * qJD(3) + t9 * qJD(4) - t16 * qJD(5) + t227; 0, 0, 0, 0, 0, 0, 0, -t241, -t187 * t293, t170, t158, 0, 0, 0, t175 * t245 - t186 * t241, t175 * t181 + t168, t37, t21, 0, 0, 0, t54 * qJD(4) - t148 * t251 + t128, t55 * qJD(4) - t146 * t251 - t129; 0, 0, 0, 0, 0, 0, 0, t213 - t261, t187 * t218 - t260, t170, t158, 0, 0, 0, t98 * qJD(4) + t186 * t213 - t266, t99 * qJD(4) + t168 - t201, t37, t21, 0, 0, 0, t25 * qJD(4) + t34 * qJD(5) + t128 - t203, t26 * qJD(4) + t33 * qJD(5) - t129 - t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t323, t181, -t245, 0, t98 * qJD(3) - t174 * t181 - t205, t99 * qJD(3) + t174 * t245 - t204, t22, t324, -t74, t75, 0, t25 * qJD(3) + t211 + t320, t26 * qJD(3) + t210 + t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t324, -t74, t75, 0, t34 * qJD(3) + t206 + t320, t33 * qJD(3) + t207 + t319; 0, 0, 0, 0, 0, 0, 0, t95 * qJD(2) + t248, t97 * qJD(2) - t249, t170, t158, 0, 0, 0, t52 * qJD(2) - t46 * qJD(4) + t186 * t248, -t51 * qJD(2) - t47 * qJD(4) - t229, t37, t21, 0, 0, 0, -t28 * qJD(2) - t3 * qJD(4) - t11 * qJD(5) - t231, t30 * qJD(2) + t5 * qJD(4) - t12 * qJD(5) + t230; 0, 0, 0, 0, 0, 0, 0, t240 + t261, t187 * t294 + t260, t170, t158, 0, 0, 0, -t94 * qJD(4) + t186 * t240 + t266, -t96 * qJD(4) + t201, t37, t21, 0, 0, 0, -t23 * qJD(4) - t31 * qJD(5) + t203, -t24 * qJD(4) - t32 * qJD(5) + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t158, 0, 0, 0, -pkin(3) * t245, -pkin(3) * t181, t37, t21, 0, 0, 0, t61 * qJD(4) - t148 * t250, t62 * qJD(4) - t146 * t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t323, t181, -t245, 0, -pkin(8) * t181 - t196, pkin(8) * t245 - t195, t22, t324, -t74, t75, 0, -t198 + t322, t197 + t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t324, -t74, t75, 0, -t194 + t322, -t193 + t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t323, 0, 0, 0, t42 * qJD(2) + t46 * qJD(3) - t233, t43 * qJD(2) + t47 * qJD(3) - t232, -t22, -t324, 0, 0, 0, t7 * qJD(2) + t3 * qJD(3) - t268, -t9 * qJD(2) - t5 * qJD(3) - t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t323, 0, 0, 0, t94 * qJD(3) + t205, t96 * qJD(3) + t204, -t22, -t324, 0, 0, 0, t23 * qJD(3) - t211, t24 * qJD(3) - t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t323, 0, 0, 0, t196, t195, -t22, -t324, 0, 0, 0, t198, -t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t304, -pkin(4) * t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242 * t304, (-t222 - t221) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t324, 0, 0, 0, t15 * qJD(2) + t11 * qJD(3) + t234, t16 * qJD(2) + t12 * qJD(3) + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t324, 0, 0, 0, t31 * qJD(3) - t206, t32 * qJD(3) - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t324, 0, 0, 0, t194, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t304, pkin(4) * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;