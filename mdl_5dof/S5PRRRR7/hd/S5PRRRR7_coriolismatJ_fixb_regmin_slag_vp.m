% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:13:01
% EndTime: 2019-12-05 17:13:10
% DurationCPUTime: 3.03s
% Computational Cost: add. (2218->159), mult. (4920->235), div. (0->0), fcn. (5364->8), ass. (0->154)
t345 = qJD(4) + qJD(5);
t352 = qJD(3) + t345;
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t177 = sin(qJ(4));
t178 = sin(qJ(3));
t181 = cos(qJ(3));
t294 = cos(qJ(4));
t151 = -t177 * t181 - t294 * t178;
t298 = pkin(6) + pkin(7);
t160 = t298 * t178;
t161 = t298 * t181;
t76 = t294 * t160 + t177 * t161;
t314 = t151 * pkin(8) - t76;
t307 = -t177 * t178 + t294 * t181;
t153 = t294 * t161;
t270 = t177 * t160;
t311 = -t153 + t270;
t71 = -pkin(8) * t307 + t311;
t349 = t352 * (-t176 * t71 - t180 * t314);
t348 = t352 * (-t176 * t314 + t180 * t71);
t236 = qJD(3) + qJD(4);
t325 = qJD(5) + t236;
t295 = cos(qJ(2));
t203 = t295 * t294;
t230 = t177 * t295;
t205 = t178 * t230;
t186 = t181 * t203 - t205;
t204 = t181 * t230;
t187 = -t178 * t203 - t204;
t297 = -t176 / 0.2e1;
t182 = t186 * t297 + t180 * t187 / 0.2e1;
t47 = t180 * t151 - t176 * t307;
t328 = t295 * t47;
t331 = -t328 / 0.2e1 + t182;
t260 = t331 * qJD(1);
t18 = t331 * qJD(2);
t179 = sin(qJ(2));
t128 = t151 * t179;
t129 = t307 * t179;
t330 = t328 / 0.2e1 + t182;
t1 = qJD(2) * t330 + t325 * (-t176 * t128 - t180 * t129);
t242 = t47 * qJD(5);
t20 = t236 * t47 + t242;
t296 = -t180 / 0.2e1;
t183 = t186 * t296 + t187 * t297;
t211 = t176 * t151 + t180 * t307;
t315 = t295 * t211;
t322 = t315 / 0.2e1 + t183;
t259 = t322 * qJD(1);
t19 = t322 * qJD(2);
t324 = t211 ^ 2 - t47 ^ 2;
t332 = t324 * qJD(2);
t323 = -t315 / 0.2e1 + t183;
t2 = qJD(2) * t323 + t325 * (-t180 * t128 + t176 * t129);
t316 = t211 * qJD(2) * t47;
t243 = qJD(5) * t211;
t21 = t236 * t211 + t243;
t173 = -t181 * pkin(3) - pkin(2);
t247 = qJD(2) * t173;
t194 = -t203 / 0.2e1;
t184 = -t204 / 0.2e1 + t178 * t194;
t224 = t295 * t151;
t96 = -t224 / 0.2e1 + t184;
t251 = t96 * qJD(1);
t308 = t151 * t247 + t251;
t290 = t178 * pkin(3);
t88 = -t173 * t151 - t290 * t307;
t306 = -t88 * qJD(2) + t251;
t291 = t151 * pkin(4);
t124 = -pkin(4) * t307 + t173;
t45 = t124 * t211;
t25 = -t291 * t47 - t45;
t305 = -t25 * qJD(2) - t259;
t44 = t124 * t47;
t24 = -t211 * t291 + t44;
t304 = -t24 * qJD(2) - t260;
t125 = t290 - t291;
t23 = -t125 * t47 + t45;
t303 = t23 * qJD(2) - t259;
t22 = -t125 * t211 - t44;
t302 = t22 * qJD(2) - t260;
t185 = t205 / 0.2e1 + t181 * t194;
t225 = t295 * t307;
t97 = t225 / 0.2e1 + t185;
t250 = t97 * qJD(1);
t301 = t236 * t76 - t250;
t202 = -t153 / 0.2e1;
t293 = pkin(3) * t177;
t285 = pkin(4) * qJD(4);
t284 = pkin(4) * qJD(5);
t235 = t294 * pkin(3);
t172 = t235 + pkin(4);
t271 = t176 * t172;
t268 = t177 * t180;
t266 = t180 * t172;
t70 = -t151 ^ 2 + t307 ^ 2;
t256 = t70 * qJD(2);
t89 = -t151 * t290 + t173 * t307;
t252 = t89 * qJD(2);
t248 = qJD(2) * t124;
t246 = qJD(2) * t181;
t245 = qJD(4) * t173;
t113 = t202 + t153 / 0.2e1;
t241 = t113 * qJD(2);
t167 = -t178 ^ 2 + t181 ^ 2;
t240 = t167 * qJD(2);
t239 = t178 * qJD(3);
t238 = t179 * qJD(2);
t237 = t181 * qJD(3);
t234 = pkin(2) * t178 * qJD(2);
t233 = pkin(2) * t246;
t223 = t294 * t176;
t219 = t307 * t247;
t217 = t178 * t246;
t216 = qJD(3) * t295;
t215 = t295 * qJD(2);
t214 = t294 * qJD(3);
t213 = t294 * qJD(4);
t212 = pkin(4) * t345;
t108 = t236 * t151;
t206 = -t235 / 0.2e1;
t188 = t206 + pkin(4) / 0.2e1 + t172 / 0.2e1;
t114 = t188 * t176;
t200 = t114 * qJD(3);
t115 = t188 * t180;
t199 = t115 * qJD(3);
t168 = t176 * t293;
t136 = t168 - t266;
t198 = t136 * qJD(3);
t137 = pkin(3) * t268 + t271;
t197 = t137 * qJD(3);
t142 = (t223 + t268) * pkin(3);
t196 = t142 * qJD(3);
t143 = t180 * t235 - t168;
t195 = t143 * qJD(3);
t190 = t211 * t248 - t259;
t189 = -t248 * t47 - t260;
t139 = t143 * qJD(4);
t138 = t142 * qJD(4);
t127 = t137 * qJD(5);
t126 = t136 * qJD(5);
t111 = t151 * t307 * qJD(2);
t110 = pkin(4) * t296 + t168 - t266 / 0.2e1 + t180 * t206;
t109 = pkin(4) * t297 - t271 / 0.2e1 + (-t268 - t223 / 0.2e1) * pkin(3);
t107 = t236 * t307;
t99 = t224 / 0.2e1 + t184;
t98 = -t225 / 0.2e1 + t185;
t80 = t97 * qJD(2);
t78 = t96 * qJD(2);
t77 = 0.2e1 * t202 + t270;
t40 = t99 * qJD(2) - t236 * t129;
t39 = t98 * qJD(2) - t236 * t128;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t238, -t215, 0, 0, 0, 0, 0, -t178 * t216 - t181 * t238, t178 * t238 - t181 * t216, 0, 0, 0, 0, 0, t236 * t99 - t238 * t307, -t151 * t238 + t236 * t98, 0, 0, 0, 0, 0, -t211 * t238 + t325 * t330, -t238 * t47 + t325 * t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178 * t215 - t179 * t237, t179 * t239 - t181 * t215, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236 * t96, -t236 * t97, 0, 0, 0, 0, 0, -t325 * t331, -t325 * t322; 0, 0, 0, 0, t178 * t237, t167 * qJD(3), 0, 0, 0, -pkin(2) * t239, -pkin(2) * t237, -t307 * t108, t236 * t70, 0, 0, 0, t88 * qJD(3) - t151 * t245, t89 * qJD(3) + t245 * t307, -t21 * t47, t325 * t324, 0, 0, 0, t22 * qJD(3) - t24 * qJD(4) - t124 * t242, t23 * qJD(3) - t25 * qJD(4) + t124 * t243; 0, 0, 0, 0, t217, t240, t237, -t239, 0, -pkin(6) * t237 - t234, pkin(6) * t239 - t233, -t111, t256, t107, t108, 0, qJD(3) * t311 + t77 * qJD(4) - t306, t252 + t301, -t316, t332, t21, t20, 0, t302 + t348, t303 + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t256, t107, t108, 0, t77 * qJD(3) + qJD(4) * t311 - t308, t219 + t301, -t316, t332, t21, t20, 0, t304 + t348, t305 + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t316, t332, t21, t20, 0, t189 + t348, t190 + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t80, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, -t217, -t240, 0, 0, 0, t234, t233, t111, -t256, 0, 0, 0, t113 * qJD(4) + t306, t250 - t252, t316, -t332, 0, 0, 0, -t302, -t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t293, -pkin(3) * t213, 0, 0, 0, 0, 0, -t138 - t127, -t139 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236 * t293 + t241, (-t214 - t213) * pkin(3), 0, 0, 0, 0, 0, t109 * qJD(5) - t138 - t196, t110 * qJD(5) - t139 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * qJD(4) - t127 - t197, t110 * qJD(4) + t126 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t80, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t256, 0, 0, 0, -t113 * qJD(3) + t308, t250 - t219, t316, -t332, 0, 0, 0, -t304, -t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t293 - t241, pkin(3) * t214, 0, 0, 0, 0, 0, -t114 * qJD(5) + t196, -t115 * qJD(5) + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176 * t284, -t180 * t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176 * t212 - t200, -t180 * t212 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, -t332, 0, 0, 0, -t189, -t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * qJD(4) + t197, t115 * qJD(4) - t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 * t285 + t200, t180 * t285 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
