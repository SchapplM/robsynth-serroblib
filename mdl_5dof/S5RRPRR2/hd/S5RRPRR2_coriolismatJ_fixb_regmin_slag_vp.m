% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR2
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
% cmat_reg [(5*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:29
% EndTime: 2019-12-05 18:28:39
% DurationCPUTime: 4.14s
% Computational Cost: add. (5399->183), mult. (10483->249), div. (0->0), fcn. (12560->8), ass. (0->169)
t339 = qJD(4) + qJD(5);
t399 = qJD(2) + t339;
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t170 = sin(pkin(9));
t171 = cos(pkin(9));
t174 = sin(qJ(2));
t177 = cos(qJ(2));
t154 = t170 * t177 + t171 * t174;
t176 = cos(qJ(4));
t152 = t176 * t154;
t153 = -t170 * t174 + t171 * t177;
t173 = sin(qJ(4));
t271 = t173 * t153;
t327 = t152 + t271;
t306 = -qJ(3) - pkin(6);
t161 = t306 * t174;
t162 = t306 * t177;
t326 = t171 * t161 + t170 * t162;
t347 = -t154 * pkin(7) + t326;
t362 = t176 * t347;
t130 = t170 * t161 - t171 * t162;
t100 = -t153 * pkin(7) - t130;
t365 = t173 * t100;
t369 = -t362 - t365;
t375 = -pkin(8) * t327 - t369;
t189 = t176 * t153 - t173 * t154;
t363 = t176 * t100;
t364 = t173 * t347;
t368 = t363 - t364;
t50 = -pkin(8) * t189 + t368;
t393 = t399 * (-t172 * t50 - t175 * t375);
t398 = t399 * (-t172 * t375 + t175 * t50);
t351 = t175 * t327;
t114 = -t351 / 0.2e1;
t207 = t351 / 0.2e1;
t370 = t114 + t207;
t394 = qJD(3) * t370 + t398;
t379 = -t172 * t327 + t175 * t189;
t241 = t379 * qJD(2);
t244 = t379 * qJD(5);
t355 = qJD(4) * t379 + t241 + t244;
t322 = qJD(3) * t379;
t382 = t379 * qJD(1);
t331 = t172 * t189;
t358 = t351 + t331;
t373 = -t358 ^ 2 + t379 ^ 2;
t378 = t373 * qJD(1);
t214 = -t177 * pkin(2) - pkin(1);
t139 = -t153 * pkin(3) + t214;
t76 = -pkin(4) * t189 + t139;
t376 = t76 * t379;
t374 = t370 * qJD(1);
t310 = -t362 / 0.2e1;
t311 = t363 / 0.2e1;
t371 = t358 * t76;
t245 = t358 * qJD(1);
t336 = t189 ^ 2 - t327 ^ 2;
t356 = t336 * qJD(1);
t308 = pkin(4) * t327;
t353 = t379 * t245;
t316 = t189 * qJD(1);
t350 = t327 * t316;
t345 = qJD(2) * t370;
t325 = 0.2e1 * t207 + t331;
t344 = qJD(2) * t325;
t343 = qJD(3) * t325;
t340 = t325 * qJD(1);
t231 = t189 * qJD(4);
t65 = t189 * qJD(2) + t231;
t321 = qJD(3) * t189;
t314 = -pkin(4) / 0.2e1;
t206 = t152 / 0.2e1;
t309 = pkin(2) * t170;
t169 = t174 * pkin(2);
t301 = pkin(4) * qJD(4);
t300 = pkin(4) * qJD(5);
t299 = qJD(2) * pkin(2);
t287 = qJD(1) * t76;
t141 = t154 * pkin(3) + t169;
t77 = t141 + t308;
t17 = -t379 * t77 + t371;
t284 = t17 * qJD(1);
t213 = t171 * pkin(2) + pkin(3);
t150 = t173 * t309 - t176 * t213;
t146 = pkin(4) - t150;
t278 = t172 * t146;
t277 = t172 * t150;
t151 = t173 * t213 + t176 * t309;
t276 = t172 * t151;
t269 = t175 * t146;
t268 = t175 * t150;
t140 = t175 * t151;
t18 = t358 * t77 + t376;
t265 = t18 * qJD(1);
t19 = t308 * t379 - t371;
t264 = t19 * qJD(1);
t20 = -t308 * t358 - t376;
t263 = t20 * qJD(1);
t55 = t214 * t169;
t251 = t55 * qJD(1);
t58 = t139 * t327 - t141 * t189;
t250 = t58 * qJD(1);
t59 = t139 * t189 + t141 * t327;
t249 = t59 * qJD(1);
t66 = t130 * t153 - t154 * t326;
t247 = t66 * qJD(1);
t240 = t358 * qJD(5);
t79 = 0.2e1 * t206 + t271;
t238 = t79 * qJD(1);
t179 = (t170 * t153 / 0.2e1 - t171 * t154 / 0.2e1) * pkin(2);
t97 = -t169 / 0.2e1 + t179;
t237 = t97 * qJD(1);
t236 = qJD(1) * t139;
t235 = qJD(1) * t177;
t121 = t206 - t152 / 0.2e1;
t234 = t121 * qJD(1);
t233 = t121 * qJD(4);
t232 = t327 * qJD(1);
t228 = t327 * qJD(4);
t126 = t153 ^ 2 + t154 ^ 2;
t227 = t126 * qJD(1);
t164 = -t174 ^ 2 + t177 ^ 2;
t226 = t164 * qJD(1);
t225 = t174 * qJD(2);
t224 = t177 * qJD(2);
t223 = pkin(1) * t174 * qJD(1);
t222 = pkin(1) * t235;
t219 = t358 * t382;
t218 = t379 * t287;
t217 = t358 * t287;
t210 = t189 * t236;
t209 = t327 * t236;
t208 = t174 * t235;
t202 = pkin(4) * t339;
t201 = t150 / 0.2e1 + pkin(4) / 0.2e1 + t146 / 0.2e1;
t98 = t140 - t277;
t199 = t98 * qJD(2);
t61 = t201 * t172;
t198 = t61 * qJD(2);
t95 = t140 + t278;
t197 = t95 * qJD(2);
t63 = t201 * t175;
t196 = t63 * qJD(2);
t94 = -t269 + t276;
t195 = t94 * qJD(2);
t99 = -t268 - t276;
t194 = t99 * qJD(2);
t192 = t339 * t370;
t23 = t310 + t362 / 0.2e1;
t188 = -t23 * qJD(1) - t150 * qJD(2);
t22 = t311 - t363 / 0.2e1;
t187 = t22 * qJD(1) + t151 * qJD(2);
t186 = qJD(2) * t327 + t79 * qJD(4);
t178 = qJD(2) * t358 + t339 * t325;
t143 = t151 * qJD(4);
t142 = t150 * qJD(4);
t96 = t169 / 0.2e1 + t179;
t93 = t99 * qJD(4);
t92 = t98 * qJD(4);
t91 = t95 * qJD(5);
t90 = t94 * qJD(5);
t64 = t175 * t314 + t276 - t269 / 0.2e1 + t268 / 0.2e1;
t62 = t172 * t314 - t140 - t278 / 0.2e1 + t277 / 0.2e1;
t57 = -t331 + 0.2e1 * t114;
t25 = 0.2e1 * t311 - t364;
t24 = -t365 + 0.2e1 * t310;
t13 = t339 * t379 + t241;
t1 = [0, 0, 0, t174 * t224, t164 * qJD(2), 0, 0, 0, -pkin(1) * t225, -pkin(1) * t224, t126 * qJD(3), t55 * qJD(2) + t66 * qJD(3), t65 * t327, (qJD(2) + qJD(4)) * t336, 0, 0, 0, t58 * qJD(2) + t139 * t228, t59 * qJD(2) + t139 * t231, t355 * t358, t399 * t373, 0, 0, 0, t17 * qJD(2) - t19 * qJD(4) + t240 * t76, t18 * qJD(2) - t20 * qJD(4) + t244 * t76; 0, 0, 0, t208, t226, t224, -t225, 0, -pkin(6) * t224 - t223, pkin(6) * t225 - t222, (-t153 * t171 - t154 * t170) * t299, t251 + t96 * qJD(3) + (-t130 * t171 + t170 * t326) * t299, t350, t356, t65, -t186, 0, qJD(2) * t368 + t25 * qJD(4) + t250, qJD(2) * t369 + t24 * qJD(4) + t249, t219, t378, t13, -t178, 0, t284 + t398, t265 + t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t96 * qJD(2) + t247, 0, 0, 0, 0, 0, t233, 0, 0, 0, 0, 0, 0, t192, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t350, t356, t65, -t79 * qJD(2) - t228, 0, t25 * qJD(2) + t121 * qJD(3) + qJD(4) * t368 + t209, t24 * qJD(2) + qJD(4) * t369 + t210, t353, t378, t355, -qJD(4) * t358 + t57 * qJD(5) - t344, 0, -t264 + t394, -t263 + t393; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, t378, t355, t57 * qJD(4) - t240 - t344, 0, t217 + t394, t218 + t393; 0, 0, 0, -t208, -t226, 0, 0, 0, t223, t222, 0, t97 * qJD(3) - t251, -t350, -t356, 0, -t233, 0, -qJD(3) * t327 - t22 * qJD(4) - t250, t23 * qJD(4) - t249 - t321, -t219, -t378, 0, -t192, 0, -qJD(3) * t358 - t284, -t265 - t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t142, 0, 0, 0, 0, 0, -t92 - t91, -t93 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, 0, 0, 0, 0, 0, -t232, -t316, 0, 0, 0, 0, 0, -t245, -t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, 0, -t143 - t187, t142 - t188, 0, 0, 0, -t374, 0, t62 * qJD(5) - t199 - t92, t64 * qJD(5) - t194 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t374, 0, t62 * qJD(4) - t197 - t91, t64 * qJD(4) + t195 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, -t97 * qJD(2) - t247, 0, 0, 0, 0, 0, t186, t65, 0, 0, 0, 0, 0, t178, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, 0, 0, 0, 0, 0, t232, t316, 0, 0, 0, 0, 0, t245, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t316, 0, 0, 0, 0, 0, t340, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, -t356, 0, t121 * qJD(2), 0, t22 * qJD(2) - t79 * qJD(3) - t209, -t23 * qJD(2) - t210 - t321, -t353, -t378, 0, qJD(5) * t370 + t345, 0, t264 - t343, t263 - t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, 0, t187, t188, 0, 0, 0, t374, 0, -t61 * qJD(5) + t199, -t63 * qJD(5) + t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, -t316, 0, 0, 0, 0, 0, -t340, -t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172 * t300, -t175 * t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t374, 0, -t172 * t202 - t198, -t175 * t202 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, -t378, 0, -qJD(4) * t370 + t345, 0, -t217 - t343, -t218 - t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t374, 0, t61 * qJD(4) + t197, t63 * qJD(4) - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, -t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t374, 0, t172 * t301 + t198, t175 * t301 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
