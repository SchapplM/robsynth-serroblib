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
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 21:23:40
% EndTime: 2021-01-15 21:23:54
% DurationCPUTime: 4.61s
% Computational Cost: add. (5441->192), mult. (10575->265), div. (0->0), fcn. (12664->8), ass. (0->175)
t347 = qJD(4) + qJD(5);
t405 = qJD(2) + t347;
t174 = sin(qJ(5));
t177 = cos(qJ(5));
t173 = sin(pkin(9));
t176 = sin(qJ(2));
t179 = cos(qJ(2));
t296 = cos(pkin(9));
t157 = t173 * t179 + t296 * t176;
t178 = cos(qJ(4));
t154 = t178 * t157;
t155 = t173 * t176 - t296 * t179;
t175 = sin(qJ(4));
t279 = t175 * t155;
t333 = t154 - t279;
t314 = -qJ(3) - pkin(6);
t163 = t314 * t176;
t164 = t314 * t179;
t206 = t296 * t163 + t173 * t164;
t187 = -t157 * pkin(7) + t206;
t335 = t178 * t187;
t132 = t173 * t163 - t296 * t164;
t100 = t155 * pkin(7) - t132;
t371 = t175 * t100;
t375 = -t335 - t371;
t381 = -pkin(8) * t333 - t375;
t192 = -t178 * t155 - t175 * t157;
t338 = t175 * t187;
t370 = t178 * t100;
t374 = t370 - t338;
t50 = -pkin(8) * t192 + t374;
t399 = t405 * (-t174 * t50 - t177 * t381);
t404 = t405 * (-t174 * t381 + t177 * t50);
t359 = t177 * t333;
t114 = -t359 / 0.2e1;
t213 = t359 / 0.2e1;
t376 = t114 + t213;
t400 = qJD(3) * t376 + t404;
t385 = -t174 * t333 + t177 * t192;
t249 = t385 * qJD(2);
t252 = t385 * qJD(5);
t363 = qJD(4) * t385 + t249 + t252;
t329 = qJD(3) * t385;
t388 = t385 * qJD(1);
t339 = t174 * t192;
t366 = t359 + t339;
t379 = -t366 ^ 2 + t385 ^ 2;
t384 = t379 * qJD(1);
t168 = -t179 * pkin(2) - pkin(1);
t141 = t155 * pkin(3) + t168;
t76 = -pkin(4) * t192 + t141;
t382 = t76 * t385;
t380 = t376 * qJD(1);
t318 = t370 / 0.2e1;
t377 = t366 * t76;
t253 = t366 * qJD(1);
t344 = t192 ^ 2 - t333 ^ 2;
t364 = t344 * qJD(1);
t317 = -t335 / 0.2e1;
t315 = pkin(4) * t333;
t361 = t385 * t253;
t323 = t192 * qJD(1);
t358 = t333 * t323;
t353 = qJD(2) * t376;
t332 = 0.2e1 * t213 + t339;
t352 = qJD(2) * t332;
t351 = qJD(3) * t332;
t348 = t332 * qJD(1);
t237 = t192 * qJD(4);
t65 = t192 * qJD(2) + t237;
t328 = qJD(3) * t192;
t321 = -pkin(4) / 0.2e1;
t212 = t154 / 0.2e1;
t316 = pkin(2) * t173;
t172 = t176 * pkin(2);
t309 = pkin(4) * qJD(4);
t308 = pkin(4) * qJD(5);
t307 = qJD(2) * pkin(2);
t294 = qJD(1) * t76;
t143 = t157 * pkin(3) + t172;
t77 = t143 + t315;
t17 = -t385 * t77 + t377;
t291 = t17 * qJD(1);
t205 = t296 * pkin(2) + pkin(3);
t152 = t175 * t316 - t178 * t205;
t148 = pkin(4) - t152;
t286 = t174 * t148;
t285 = t174 * t152;
t153 = t175 * t205 + t178 * t316;
t284 = t174 * t153;
t277 = t177 * t148;
t276 = t177 * t152;
t142 = t177 * t153;
t18 = t366 * t77 + t382;
t273 = t18 * qJD(1);
t19 = t315 * t385 - t377;
t272 = t19 * qJD(1);
t20 = -t315 * t366 - t382;
t271 = t20 * qJD(1);
t55 = t168 * t172;
t259 = t55 * qJD(1);
t58 = t141 * t333 - t143 * t192;
t258 = t58 * qJD(1);
t59 = t141 * t192 + t143 * t333;
t257 = t59 * qJD(1);
t66 = -t132 * t155 - t157 * t206;
t255 = t66 * qJD(1);
t248 = t366 * qJD(5);
t79 = 0.2e1 * t212 - t279;
t246 = t79 * qJD(1);
t180 = (-t173 * t155 / 0.2e1 - t296 * t157 / 0.2e1) * pkin(2);
t97 = -t172 / 0.2e1 + t180;
t245 = t97 * qJD(1);
t244 = qJD(1) * t141;
t243 = qJD(1) * t179;
t116 = t155 * t172 + t168 * t157;
t242 = t116 * qJD(1);
t117 = -t168 * t155 + t157 * t172;
t241 = t117 * qJD(1);
t123 = t212 - t154 / 0.2e1;
t240 = t123 * qJD(1);
t239 = t123 * qJD(4);
t238 = t333 * qJD(1);
t234 = t333 * qJD(4);
t128 = t155 ^ 2 + t157 ^ 2;
t233 = t128 * qJD(1);
t232 = t155 * qJD(1);
t231 = t157 * qJD(1);
t166 = -t176 ^ 2 + t179 ^ 2;
t230 = t166 * qJD(1);
t229 = t176 * qJD(2);
t228 = t179 * qJD(2);
t227 = pkin(1) * t176 * qJD(1);
t226 = pkin(1) * t243;
t223 = t366 * t388;
t222 = t385 * t294;
t221 = t366 * t294;
t216 = t192 * t244;
t215 = t333 * t244;
t214 = t176 * t243;
t207 = pkin(4) * t347;
t204 = t152 / 0.2e1 + pkin(4) / 0.2e1 + t148 / 0.2e1;
t98 = t142 - t285;
t202 = t98 * qJD(2);
t61 = t204 * t174;
t201 = t61 * qJD(2);
t95 = t142 + t286;
t200 = t95 * qJD(2);
t63 = t204 * t177;
t199 = t63 * qJD(2);
t94 = -t277 + t284;
t198 = t94 * qJD(2);
t99 = -t276 - t284;
t197 = t99 * qJD(2);
t195 = t347 * t376;
t23 = t317 + t335 / 0.2e1;
t191 = -t23 * qJD(1) - t152 * qJD(2);
t22 = t318 - t370 / 0.2e1;
t190 = t22 * qJD(1) + t153 * qJD(2);
t189 = qJD(2) * t333 + t79 * qJD(4);
t181 = qJD(2) * t366 + t332 * t347;
t145 = t153 * qJD(4);
t144 = t152 * qJD(4);
t96 = t172 / 0.2e1 + t180;
t93 = t99 * qJD(4);
t92 = t98 * qJD(4);
t91 = t95 * qJD(5);
t90 = t94 * qJD(5);
t64 = t177 * t321 + t284 - t277 / 0.2e1 + t276 / 0.2e1;
t62 = t174 * t321 - t142 - t286 / 0.2e1 + t285 / 0.2e1;
t57 = -t339 + 0.2e1 * t114;
t25 = 0.2e1 * t318 - t338;
t24 = -t371 + 0.2e1 * t317;
t13 = t347 * t385 + t249;
t1 = [0, 0, 0, t176 * t228, t166 * qJD(2), 0, 0, 0, -pkin(1) * t229, -pkin(1) * t228, t116 * qJD(2), t117 * qJD(2), t128 * qJD(3), t55 * qJD(2) + t66 * qJD(3), t65 * t333, (qJD(2) + qJD(4)) * t344, 0, 0, 0, t58 * qJD(2) + t141 * t234, t59 * qJD(2) + t141 * t237, t363 * t366, t405 * t379, 0, 0, 0, t17 * qJD(2) - t19 * qJD(4) + t248 * t76, t18 * qJD(2) - t20 * qJD(4) + t252 * t76; 0, 0, 0, t214, t230, t228, -t229, 0, -pkin(6) * t228 - t227, pkin(6) * t229 - t226, -qJD(2) * t132 + t242, -qJD(2) * t206 + t241, (t296 * t155 - t157 * t173) * t307, t259 + (-t132 * t296 + t173 * t206) * t307 + t96 * qJD(3), t358, t364, t65, -t189, 0, qJD(2) * t374 + t25 * qJD(4) + t258, qJD(2) * t375 + t24 * qJD(4) + t257, t223, t384, t13, -t181, 0, t291 + t404, t273 + t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t96 * qJD(2) + t255, 0, 0, 0, 0, 0, t239, 0, 0, 0, 0, 0, 0, t195, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t358, t364, t65, -t79 * qJD(2) - t234, 0, t25 * qJD(2) + t123 * qJD(3) + qJD(4) * t374 + t215, t24 * qJD(2) + qJD(4) * t375 + t216, t361, t384, t363, -qJD(4) * t366 + t57 * qJD(5) - t352, 0, -t272 + t400, -t271 + t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t361, t384, t363, t57 * qJD(4) - t248 - t352, 0, t221 + t400, t222 + t399; 0, 0, 0, -t214, -t230, 0, 0, 0, t227, t226, -t157 * qJD(3) - t242, t155 * qJD(3) - t241, 0, t97 * qJD(3) - t259, -t358, -t364, 0, -t239, 0, -qJD(3) * t333 - t22 * qJD(4) - t258, t23 * qJD(4) - t257 - t328, -t223, -t384, 0, -t195, 0, -qJD(3) * t366 - t291, -t273 - t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, t144, 0, 0, 0, 0, 0, -t92 - t91, -t93 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, t232, 0, t245, 0, 0, 0, 0, 0, -t238, -t323, 0, 0, 0, 0, 0, -t253, -t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, 0, -t145 - t190, t144 - t191, 0, 0, 0, -t380, 0, t62 * qJD(5) - t202 - t92, t64 * qJD(5) - t197 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, 0, t62 * qJD(4) - t200 - t91, t64 * qJD(4) + t198 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157 * qJD(2), -t155 * qJD(2), -t233, -t97 * qJD(2) - t255, 0, 0, 0, 0, 0, t189, t65, 0, 0, 0, 0, 0, t181, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, -t232, 0, -t245, 0, 0, 0, 0, 0, t238, t323, 0, 0, 0, 0, 0, t253, t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, t323, 0, 0, 0, 0, 0, t348, t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t358, -t364, 0, t123 * qJD(2), 0, t22 * qJD(2) - t79 * qJD(3) - t215, -t23 * qJD(2) - t216 - t328, -t361, -t384, 0, qJD(5) * t376 + t353, 0, t272 - t351, t271 - t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, t190, t191, 0, 0, 0, t380, 0, -t61 * qJD(5) + t202, -t63 * qJD(5) + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, -t323, 0, 0, 0, 0, 0, -t348, -t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174 * t308, -t177 * t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t380, 0, -t174 * t207 - t201, -t177 * t207 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361, -t384, 0, -qJD(4) * t376 + t353, 0, -t221 - t351, -t222 - t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t380, 0, t61 * qJD(4) + t200, t63 * qJD(4) - t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t348, -t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, 0, t174 * t309 + t201, t177 * t309 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
