% Calculate inertial parameters regressor of coriolis matrix for
% S5RRPRP3
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
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:18
% EndTime: 2019-12-31 19:51:23
% DurationCPUTime: 2.49s
% Computational Cost: add. (3697->235), mult. (7250->254), div. (0->0), fcn. (7313->6), ass. (0->175)
t193 = sin(pkin(8));
t194 = cos(pkin(8));
t195 = sin(qJ(4));
t296 = cos(qJ(4));
t166 = t296 * t193 + t195 * t194;
t197 = cos(qJ(2));
t294 = t197 * pkin(1);
t136 = t166 * t294;
t204 = -t195 * t193 + t296 * t194;
t137 = t204 * t294;
t210 = t136 * t166 + t137 * t204;
t302 = t210 * qJD(1);
t160 = t204 ^ 2;
t299 = t166 ^ 2;
t309 = t160 + t299;
t320 = t309 * qJD(3);
t325 = qJD(2) * t210 + t320;
t329 = t302 + t325;
t245 = qJD(1) + qJD(2);
t97 = t299 - t160;
t328 = t245 * t97;
t327 = t97 * qJD(4);
t326 = t320 - t302;
t311 = t245 * t166;
t324 = t204 * t311;
t323 = t245 * t309;
t190 = t194 * pkin(7);
t196 = sin(qJ(2));
t295 = t196 * pkin(1);
t234 = qJ(3) + t295;
t159 = t194 * t234 + t190;
t208 = (-pkin(7) - t234) * t193;
t103 = t195 * t159 - t296 * t208;
t104 = t296 * t159 + t195 * t208;
t180 = t194 * qJ(3) + t190;
t235 = (-pkin(7) - qJ(3)) * t193;
t126 = t195 * t180 - t296 * t235;
t127 = t296 * t180 + t195 * t235;
t188 = t295 / 0.2e1;
t198 = t188 + (-t126 / 0.2e1 - t103 / 0.2e1) * t166 - (t127 / 0.2e1 + t104 / 0.2e1) * t204;
t280 = t104 * t204;
t282 = t103 * t166;
t212 = t280 + t282;
t300 = t212 * qJD(1);
t322 = qJD(2) * t198 - t300;
t276 = t127 * t204;
t277 = t126 * t166;
t200 = t276 / 0.2e1 + t280 / 0.2e1 + t277 / 0.2e1 + t282 / 0.2e1 + t188;
t321 = qJD(2) * t200 + t300;
t116 = pkin(4) * t166 - qJ(5) * t204;
t319 = t245 * t116;
t187 = -t194 * pkin(3) - pkin(2);
t179 = t187 - t294;
t236 = t187 / 0.2e1 + t179 / 0.2e1;
t316 = t236 * t166;
t221 = -pkin(4) * t204 - t166 * qJ(5);
t106 = t187 + t221;
t99 = t106 - t294;
t240 = t106 / 0.2e1 + t99 / 0.2e1;
t315 = t240 * t166;
t312 = t245 * t204;
t191 = t193 ^ 2;
t192 = t194 ^ 2;
t183 = t191 + t192;
t310 = t245 * t183;
t305 = qJD(3) * t198;
t304 = qJD(3) * t212;
t303 = t200 * qJD(3);
t211 = t276 + t277;
t301 = t211 * qJD(3);
t138 = t183 * t234;
t216 = qJD(1) * t198 - qJD(2) * t211;
t293 = pkin(1) * qJD(1);
t292 = pkin(1) * qJD(2);
t9 = t99 * t116;
t289 = t9 * qJD(1);
t288 = t99 * t204;
t287 = t99 * t166;
t155 = t166 * qJD(5);
t284 = qJD(4) * t116 - t155;
t279 = t106 * t204;
t278 = t106 * t166;
t16 = t116 * t106;
t275 = t179 * t204;
t274 = t187 * t204;
t213 = t103 * t136 + t104 * t137;
t23 = t99 * t295 + t213;
t273 = t23 * qJD(1);
t25 = t179 * t295 + t213;
t272 = t25 * qJD(1);
t84 = t116 * t204;
t31 = -t84 + t287;
t269 = t31 * qJD(1);
t85 = t116 * t166;
t32 = -t85 - t288;
t268 = t32 * qJD(1);
t265 = t221 * qJD(4) + qJD(5) * t204;
t125 = t204 * t155;
t243 = t196 * t292;
t141 = t204 * t243;
t264 = -t141 + t125;
t231 = t183 * t197;
t158 = pkin(1) * t231;
t181 = t183 * qJD(3);
t263 = t158 * qJD(2) + t181;
t154 = t299 * qJD(5);
t229 = t166 * t243;
t262 = t154 - t229;
t230 = t296 * t294;
t219 = t230 / 0.2e1;
t244 = t195 * t294;
t226 = -t244 / 0.2e1;
t261 = t193 * t226 + t194 * t219;
t225 = t244 / 0.2e1;
t260 = t193 * t219 + t194 * t225;
t220 = -t230 / 0.2e1;
t259 = t193 * t220 + t194 * t226;
t258 = t193 * t225 + t194 * t220;
t257 = qJD(4) * qJ(5);
t256 = t103 * qJD(4);
t100 = t104 * qJD(4);
t105 = (-pkin(2) - t294) * t295 + t294 * t138;
t255 = t105 * qJD(1);
t254 = t126 * qJD(4);
t119 = t127 * qJD(4);
t253 = t138 * qJD(1);
t252 = t158 * qJD(1);
t251 = t204 * qJD(3);
t250 = t204 * qJD(4);
t249 = t166 * qJD(1);
t248 = t166 * qJD(2);
t247 = t166 * qJD(3);
t246 = t166 * qJD(4);
t242 = t196 * t293;
t241 = t99 * t249;
t239 = t204 * t246;
t238 = qJD(1) * t275;
t237 = t179 * t249;
t233 = t136 * t126 + t137 * t127;
t232 = pkin(1) * t245;
t178 = t183 * qJ(3);
t228 = t166 * t242;
t227 = t193 * t242;
t224 = t196 * t232;
t206 = t137 * qJ(5) / 0.2e1 - t136 * pkin(4) / 0.2e1;
t1 = -t240 * t116 + t206;
t218 = -t1 * qJD(1) + t16 * qJD(2);
t17 = t259 + t84 - t315;
t33 = -t84 + t278;
t215 = -t17 * qJD(1) + t33 * qJD(2);
t18 = t204 * t240 + t261 + t85;
t34 = -t85 - t279;
t214 = t18 * qJD(1) - t34 * qJD(2);
t199 = t178 + (t191 / 0.2e1 + t192 / 0.2e1) * t295;
t102 = -t295 / 0.2e1 + t199;
t209 = t102 * qJD(1) + t178 * qJD(2);
t205 = -t278 / 0.2e1 - t287 / 0.2e1;
t26 = t260 + t315;
t203 = t26 * qJD(1) + t106 * t248;
t48 = -t204 * t236 + t258;
t202 = t48 * qJD(1) - qJD(2) * t274;
t47 = t259 - t316;
t201 = -t47 * qJD(1) + t187 * t248;
t182 = t193 * t243;
t139 = t204 * t242;
t114 = t245 * t299;
t101 = t188 + t199;
t80 = t116 * qJD(3);
t50 = t259 + t316;
t49 = t274 / 0.2e1 + t275 / 0.2e1 + t258;
t27 = t205 + t260;
t20 = -t85 - t279 / 0.2e1 - t288 / 0.2e1 + t261;
t19 = -t84 - t205 + t259;
t2 = t16 / 0.2e1 + t9 / 0.2e1 + t206;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, -t197 * t292, 0, 0, 0, 0, 0, 0, 0, 0, -t194 * t243, t182, t263, t105 * qJD(2) + t138 * qJD(3), t239, -t327, 0, -t239, 0, 0, t179 * t246 - t141, t179 * t250 + t229, t325, t25 * qJD(2) + t304, t239, 0, t327, 0, 0, -t239, t31 * qJD(4) + t264, t325, t32 * qJD(4) + t262, t23 * qJD(2) + t9 * qJD(4) - t155 * t99 + t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, -t197 * t232, 0, 0, 0, 0, 0, 0, 0, 0, -t194 * t224, t182 + t227, t252 + t263, t255 + t101 * qJD(3) + (-pkin(2) * t196 + qJ(3) * t231) * t292, t239, -t327, 0, -t239, 0, 0, t50 * qJD(4) - t139 - t141, t49 * qJD(4) + t228 + t229, t329, t272 + (t187 * t295 + t233) * qJD(2) + t303, t239, 0, t327, 0, 0, -t239, t19 * qJD(4) - t139 + t264, t329, t20 * qJD(4) - t228 + t262, t273 + (t106 * t295 + t233) * qJD(2) + t303 + t2 * qJD(4) + t27 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t101 * qJD(2) + t253, 0, 0, 0, 0, 0, 0, 0, 0, t323, t321, 0, 0, 0, 0, 0, 0, 0, t323, 0, t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, -t328, t250, -t324, -t246, 0, t50 * qJD(2) - t100 + t237, t49 * qJD(2) + t238 + t256, 0, 0, t324, t250, t328, 0, t246, -t324, t19 * qJD(2) - t100 + t269, t265, t20 * qJD(2) - t256 + t268, t289 + t2 * qJD(2) + (-t104 * pkin(4) - t103 * qJ(5)) * qJD(4) + t104 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t250, t114, t27 * qJD(2) + t100 - t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t197 * t293, 0, 0, 0, 0, 0, 0, 0, 0, t194 * t242, -t227, t181 - t252, t102 * qJD(3) - t255, t239, -t327, 0, -t239, 0, 0, -t47 * qJD(4) + t139, -t48 * qJD(4) - t228, t326, -t272 - t305, t239, 0, t327, 0, 0, -t239, -t17 * qJD(4) + t125 + t139, t326, -t18 * qJD(4) + t154 + t228, -t1 * qJD(4) - t26 * qJD(5) - t273 - t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, t178 * qJD(3), t239, -t327, 0, -t239, 0, 0, t187 * t246, t187 * t250, t320, t301, t239, 0, t327, 0, 0, -t239, t33 * qJD(4) + t125, t320, t34 * qJD(4) + t154, t16 * qJD(4) - t106 * t155 + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t209, 0, 0, 0, 0, 0, 0, 0, 0, t323, -t216, 0, 0, 0, 0, 0, 0, 0, t323, 0, -t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, -t328, t250, -t324, -t246, 0, -t119 + t201, -t202 + t254, 0, 0, t324, t250, t328, 0, t246, -t324, -t119 + t215, t265, -t214 - t254, (-t127 * pkin(4) - t126 * qJ(5)) * qJD(4) + t127 * qJD(5) + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t250, t114, -t203 + t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, -t102 * qJD(2) - t253, 0, 0, 0, 0, 0, 0, t246, t250, -t323, t322, 0, 0, 0, 0, 0, 0, t246, -t323, -t250, t284 + t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, -t209, 0, 0, 0, 0, 0, 0, t246, t250, -t323, t216, 0, 0, 0, 0, 0, 0, t246, -t323, -t250, t216 + t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, t312, 0, 0, 0, 0, 0, 0, 0, 0, t311, 0, -t312, t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, t328, 0, t324, 0, 0, t47 * qJD(2) - t237 - t247, t48 * qJD(2) - t238 - t251, 0, 0, -t324, 0, -t328, 0, 0, t324, t17 * qJD(2) - t247 - t269, 0, t18 * qJD(2) + t251 - t268, t1 * qJD(2) - t289 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, t328, 0, t324, 0, 0, -t247 - t201, -t251 + t202, 0, 0, -t324, 0, -t328, 0, 0, t324, -t247 - t215, 0, t251 + t214, -t218 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t311, -t312, 0, 0, 0, 0, 0, 0, 0, 0, -t311, 0, t312, -t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, 0, -t114, t26 * qJD(2) + t241 + t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, 0, -t114, t247 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
