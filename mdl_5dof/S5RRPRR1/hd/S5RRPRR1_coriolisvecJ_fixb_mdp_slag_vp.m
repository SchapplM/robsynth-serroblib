% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:40
% EndTime: 2019-12-05 18:25:47
% DurationCPUTime: 2.44s
% Computational Cost: add. (1370->233), mult. (3402->334), div. (0->0), fcn. (2280->6), ass. (0->120)
t246 = qJD(2) + qJD(4);
t251 = sin(qJ(4));
t254 = cos(qJ(2));
t320 = cos(qJ(4));
t284 = qJD(1) * t320;
t252 = sin(qJ(2));
t300 = qJD(1) * t252;
t325 = -t251 * t300 + t254 * t284;
t326 = t325 * t246;
t221 = t251 * t254 + t320 * t252;
t197 = t246 * t221;
t191 = t197 * qJD(1);
t255 = pkin(1) + pkin(2);
t213 = qJD(5) - t325;
t324 = t213 ^ 2;
t247 = t252 ^ 2;
t248 = t254 ^ 2;
t323 = MDP(5) * (t247 - t248);
t319 = pkin(3) + qJ(3);
t227 = t319 * t252;
t222 = qJD(1) * t227;
t249 = qJD(2) * pkin(1);
t209 = qJD(2) * pkin(2) - t222 + t249;
t281 = qJD(2) * t319;
t210 = qJD(3) * t254 - t252 * t281;
t207 = t210 * qJD(1);
t296 = qJD(3) * t252;
t211 = -t254 * t281 - t296;
t208 = t211 * qJD(1);
t228 = t319 * t254;
t223 = qJD(1) * t228;
t295 = qJD(4) * t251;
t262 = t320 * t207 + t251 * t208 - t223 * t295;
t283 = qJD(4) * t320;
t169 = t209 * t283 + t262;
t288 = t320 * t223;
t189 = t251 * t209 + t288;
t279 = t251 * t207 - t320 * t208;
t170 = qJD(4) * t189 + t279;
t265 = -t320 * t227 - t251 * t228;
t173 = t265 * qJD(4) + t320 * t210 + t251 * t211;
t289 = t320 * t209;
t305 = t251 * t223;
t188 = -t289 + t305;
t299 = qJD(1) * t254;
t216 = -t251 * t299 - t252 * t284;
t229 = t255 * t254;
t219 = -qJD(1) * t229 + qJD(3);
t195 = pkin(4) * t216 + t219;
t264 = -t251 * t252 + t320 * t254;
t196 = t246 * t264;
t203 = -t251 * t227 + t320 * t228;
t205 = -pkin(4) * t221 - t229;
t322 = t170 * t221 + t188 * t196 - t203 * t191 - (qJD(5) * t205 + t173) * t213 + (qJD(5) * t195 + t169) * t264;
t321 = 0.2e1 * t248;
t250 = sin(qJ(5));
t294 = qJD(5) * t250;
t253 = cos(qJ(5));
t293 = qJD(5) * t253;
t302 = t246 * t293 + t253 * t326;
t175 = t216 * t294 + t302;
t318 = t175 * t250;
t317 = t188 * t221;
t316 = t326 * t250;
t308 = t216 * t250;
t199 = -t253 * t246 - t308;
t315 = t199 * t213;
t314 = t199 * t216;
t270 = t216 * t253 - t246 * t250;
t313 = t270 * t213;
t312 = t270 * t216;
t311 = t205 * t191;
t309 = t216 * t246;
t307 = t219 * t325;
t306 = t250 * t191;
t257 = qJD(1) ^ 2;
t304 = t252 * t257;
t187 = t253 * t191;
t303 = t253 * t221;
t292 = qJD(1) * qJD(2);
t282 = t252 * t292;
t239 = pkin(1) * t282;
t217 = pkin(2) * t282 + t239;
t224 = t255 * qJD(2) * t252;
t225 = t255 * t300;
t297 = qJD(2) * t254;
t291 = pkin(1) * t299;
t290 = t255 * t320;
t287 = qJ(3) * t297;
t226 = -qJ(3) * t300 + t249;
t285 = t226 * t297;
t280 = t188 * (-t213 - t325);
t278 = t253 * t213;
t234 = t251 * t255 + pkin(4);
t275 = -pkin(4) * t325 + qJD(5) * t234 + t225;
t183 = t246 * pkin(4) + t189;
t172 = t183 * t253 + t195 * t250;
t274 = t170 * t250 - t172 * t216 + t188 * t293;
t193 = -t251 * t222 + t288;
t273 = t255 * t295 - t193;
t272 = t183 * t250 - t195 * t253;
t271 = -t188 * t325 - t191 * t234;
t269 = t187 + (t250 * t325 - t294) * t213;
t194 = -t320 * t222 - t305;
t268 = -t255 * t283 + t194;
t267 = -t170 * t253 + t188 * t294 - t216 * t272;
t266 = t219 * t216 - t279;
t263 = t253 * t196 - t221 * t294;
t260 = -t324 * t253 - t306;
t176 = -qJD(5) * t270 + t316;
t259 = ((t175 - t315) * t253 + (-t176 + t313) * t250) * MDP(21) + (-t270 * t278 + t318) * MDP(20) + (t269 - t314) * MDP(23) + (t213 * t278 + t306 - t312) * MDP(22) + (-t309 - t191) * MDP(16) + (t216 ^ 2 - t325 ^ 2) * MDP(14) + (MDP(13) * t325 + t213 * MDP(24)) * t216;
t258 = qJ(3) ^ 2;
t256 = qJD(2) ^ 2;
t233 = qJD(3) - t291;
t218 = (-t287 - t296) * qJD(1);
t182 = -pkin(4) * t196 + t224;
t180 = -pkin(4) * t326 + t217;
t179 = t253 * t180;
t174 = t203 * qJD(4) + t251 * t210 - t320 * t211;
t1 = [-0.2e1 * t292 * t323 - t256 * t252 * MDP(7) + (-t285 - t218 * t252 + (-0.2e1 * t252 * t287 + (t247 + t321) * qJD(3)) * qJD(1)) * MDP(11) + ((qJD(1) * qJD(3) * t321 - t285) * qJ(3) + (-qJ(3) * t218 - qJD(3) * t226 + (-0.2e1 * t258 * t299 + (t233 - t291) * pkin(1)) * qJD(2)) * t252) * MDP(12) + (-t196 * t216 + t221 * t326) * MDP(13) + (-t191 * t221 + t196 * t325 + t197 * t216 + t264 * t326) * MDP(14) + (-t191 * t229 + t197 * t219 - t217 * t264 - t224 * t325) * MDP(18) + (t196 * t219 - t216 * t224 + t217 * t221 - t229 * t326) * MDP(19) + (t175 * t303 - t263 * t270) * MDP(20) + ((-t199 * t253 + t250 * t270) * t196 + (-t318 - t176 * t253 + (t199 * t250 + t253 * t270) * qJD(5)) * t221) * MDP(21) + (-t175 * t264 + t191 * t303 - t197 * t270 + t213 * t263) * MDP(22) + (-t221 * t306 + t176 * t264 - t197 * t199 + (-t250 * t196 - t221 * t293) * t213) * MDP(23) + (-t191 * t264 + t197 * t213) * MDP(24) + (-t272 * t197 + t174 * t199 - t265 * t176 - t179 * t264 + (t182 * t213 + t311 + (t183 * t264 - t203 * t213 + t317) * qJD(5)) * t253 + t322 * t250) * MDP(25) + (-t172 * t197 - t174 * t270 - t265 * t175 + (-(-qJD(5) * t203 + t182) * t213 - t311 + (-qJD(5) * t183 + t180) * t264 - qJD(5) * t317) * t250 + t322 * t253) * MDP(26) + (0.2e1 * MDP(4) * t282 + t256 * MDP(6)) * t254 + (t196 * MDP(15) - t197 * MDP(16) - t174 * MDP(18) - t173 * MDP(19)) * t246; t259 + (-t175 * t290 + t271 * t253 - t273 * t270 + (t275 * t250 + t268 * t253) * t213 + t274) * MDP(26) + (-t176 * t290 + t271 * t250 + t273 * t199 + (t268 * t250 - t275 * t253) * t213 + t267) * MDP(25) + (t193 * t246 + t225 * t325 + (-t288 + (-t246 * t255 - t209) * t251) * qJD(4) + t266) * MDP(18) + (t194 * t246 - t307 + t225 * t216 + (-t246 * t290 - t289) * qJD(4) - t262) * MDP(19) + (-t233 * t300 + t218) * pkin(1) * MDP(12) + t257 * t323 + ((qJ(3) * qJD(1) * t226 + t258 * t304) * MDP(12) + (qJ(3) * t304 + (t226 - t249) * qJD(1)) * MDP(11) - MDP(4) * t304) * t254; (-t247 - t248) * t257 * MDP(11) + (-qJ(3) * t248 * t257 + t226 * t300 + t239) * MDP(12) + (t191 - t309) * MDP(18) + 0.2e1 * t326 * MDP(19) + (t269 + t314) * MDP(25) + (t260 - t312) * MDP(26); (t266 + (-qJD(4) + t246) * t189) * MDP(18) + (-t188 * t246 - t169 - t307) * MDP(19) + (t260 * pkin(4) - t189 * t199 + t250 * t280 + t267) * MDP(25) + (t189 * t270 + t253 * t280 + (t324 * t250 - t187) * pkin(4) + t274) * MDP(26) + t259; -t270 * t199 * MDP(20) + (-t199 ^ 2 + t270 ^ 2) * MDP(21) + (t302 + t315) * MDP(22) + (-t313 - t316) * MDP(23) + t191 * MDP(24) + (-t169 * t250 + t172 * t213 + t188 * t270 + t179) * MDP(25) + (-t169 * t253 - t180 * t250 + t188 * t199 - t213 * t272) * MDP(26) + (MDP(22) * t308 + t270 * MDP(23) - MDP(25) * t172 + t272 * MDP(26)) * qJD(5);];
tauc = t1;
