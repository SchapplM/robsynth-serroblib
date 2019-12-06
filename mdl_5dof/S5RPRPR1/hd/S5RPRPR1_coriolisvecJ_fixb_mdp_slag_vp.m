% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:48
% EndTime: 2019-12-05 17:47:52
% DurationCPUTime: 1.24s
% Computational Cost: add. (1049->177), mult. (2329->251), div. (0->0), fcn. (1574->6), ass. (0->93)
t249 = sin(pkin(8));
t250 = cos(pkin(8));
t252 = sin(qJ(3));
t254 = cos(qJ(3));
t261 = t249 * t254 + t250 * t252;
t219 = t261 * qJD(1);
t253 = cos(qJ(5));
t209 = t253 * t219;
t260 = t249 * t252 - t250 * t254;
t268 = qJD(1) * qJD(3);
t213 = t260 * t268;
t214 = t261 * t268;
t278 = qJD(1) * t254;
t279 = qJD(1) * t252;
t222 = -t249 * t279 + t250 * t278;
t251 = sin(qJ(5));
t273 = qJD(5) * t251;
t152 = -qJD(5) * t209 + t251 * t213 - t253 * t214 - t222 * t273;
t262 = -t219 * t251 + t253 * t222;
t153 = t262 * qJD(5) - t253 * t213 - t214 * t251;
t177 = t222 * t251 + t209;
t244 = qJD(3) + qJD(5);
t288 = t177 * t244;
t289 = t262 * t244;
t303 = t177 * t262 * MDP(16) + (-t153 + t289) * MDP(19) + (-t177 ^ 2 + t262 ^ 2) * MDP(17) + (t152 + t288) * MDP(18);
t276 = qJD(3) * t254;
t269 = pkin(3) * t276 + qJD(2);
t302 = qJ(2) * (t254 * MDP(12) - t252 * MDP(13)) - t254 * t252 * MDP(7) + (t252 ^ 2 - t254 ^ 2) * MDP(8);
t255 = -pkin(1) - pkin(6);
t233 = t255 * qJD(1) + qJD(2);
t218 = -qJ(4) * t278 + t254 * t233;
t206 = qJD(3) * pkin(3) + t218;
t217 = -qJ(4) * t279 + t233 * t252;
t287 = t250 * t217;
t169 = t249 * t206 + t287;
t292 = pkin(7) * t219;
t161 = t169 - t292;
t230 = pkin(3) * t279 + qJD(1) * qJ(2) + qJD(4);
t189 = pkin(4) * t219 + t230;
t301 = t161 * t273 + t189 * t177;
t297 = MDP(12) * t252 + MDP(13) * t254;
t296 = qJD(5) - t244;
t274 = qJD(4) * t254;
t277 = qJD(3) * t252;
t190 = -t233 * t277 + (qJ(4) * t277 - t274) * qJD(1);
t275 = qJD(4) * t252;
t191 = t233 * t276 + (-qJ(4) * t276 - t275) * qJD(1);
t164 = t250 * t190 - t191 * t249;
t155 = pkin(7) * t214 + t164;
t165 = t249 * t190 + t250 * t191;
t156 = pkin(7) * t213 + t165;
t294 = t253 * t155 - t251 * t156 - t189 * t262;
t293 = pkin(3) * t249;
t291 = pkin(7) * t222;
t202 = t249 * t217;
t286 = t252 * pkin(3) + qJ(2);
t285 = qJ(4) - t255;
t215 = t285 * t277 - t274;
t232 = t285 * t254;
t216 = -qJD(3) * t232 - t275;
t171 = t249 * t215 + t250 * t216;
t175 = t250 * t218 - t202;
t231 = t285 * t252;
t188 = -t250 * t231 - t249 * t232;
t284 = t269 * qJD(1);
t266 = -qJ(2) * MDP(6) - MDP(5);
t168 = t250 * t206 - t202;
t170 = t250 * t215 - t216 * t249;
t174 = -t218 * t249 - t287;
t187 = t231 * t249 - t250 * t232;
t160 = qJD(3) * pkin(4) + t168 - t291;
t263 = -t251 * t160 - t253 * t161;
t182 = -t251 * t261 - t253 * t260;
t181 = -t251 * t260 + t253 * t261;
t220 = t249 * t277 - t250 * t276;
t221 = t261 * qJD(3);
t259 = -t164 * t260 + t165 * t261 - t168 * t221 - t169 * t220;
t157 = -t181 * qJD(5) + t220 * t251 - t253 * t221;
t158 = t182 * qJD(5) - t253 * t220 - t221 * t251;
t257 = qJD(1) ^ 2;
t256 = qJD(3) ^ 2;
t239 = pkin(3) * t250 + pkin(4);
t207 = pkin(4) * t261 + t286;
t195 = pkin(3) * t278 + pkin(4) * t222;
t192 = -pkin(4) * t220 + t269;
t183 = -pkin(4) * t213 + t284;
t173 = -pkin(7) * t261 + t188;
t172 = pkin(7) * t260 + t187;
t167 = t175 - t291;
t166 = t174 + t292;
t163 = pkin(7) * t220 + t171;
t162 = pkin(7) * t221 + t170;
t1 = [(-t170 * t222 - t171 * t219 + t187 * t214 + t188 * t213 - t259) * MDP(14) + (t164 * t187 + t165 * t188 + t168 * t170 + t169 * t171 + t230 * t269 + t284 * t286) * MDP(15) + (t152 * t182 + t157 * t262) * MDP(16) + (-t152 * t181 - t153 * t182 - t157 * t177 - t158 * t262) * MDP(17) + (t207 * t153 + t189 * t158 + t192 * t177 + t183 * t181) * MDP(21) + (t207 * t152 + t189 * t157 + t183 * t182 + t192 * t262) * MDP(22) + (t157 * MDP(18) - t158 * MDP(19) + (t162 * t253 - t163 * t251) * MDP(21) + (-t162 * t251 - t163 * t253) * MDP(22) + ((-t172 * t251 - t173 * t253) * MDP(21) + (-t172 * t253 + t173 * t251) * MDP(22)) * qJD(5)) * t244 + ((-MDP(13) * t255 - MDP(10)) * t254 + (-MDP(12) * t255 - MDP(9)) * t252) * t256 + (0.2e1 * (-t266 + t297) * qJD(2) + 0.2e1 * t302 * qJD(3)) * qJD(1); (t213 * t261 - t214 * t260 + t219 * t220 + t221 * t222) * MDP(14) + (-qJD(1) * t230 + t259) * MDP(15) + (-qJD(1) * t177 + t157 * t244) * MDP(21) + (-qJD(1) * t262 - t158 * t244) * MDP(22) + t266 * t257 + t297 * (-t256 - t257); ((t169 + t174) * t222 - (t168 - t175) * t219 + (t213 * t249 + t214 * t250) * pkin(3)) * MDP(14) + (-t168 * t174 - t169 * t175 + (t164 * t250 + t165 * t249 - t230 * t278) * pkin(3)) * MDP(15) + (-t195 * t177 - (t166 * t253 - t167 * t251) * t244 + ((-t239 * t251 - t253 * t293) * t244 + t263) * qJD(5) + t294) * MDP(21) + (-t253 * t156 - t251 * t155 - t195 * t262 + (t166 * t251 + t167 * t253) * t244 + (-(t239 * t253 - t251 * t293) * t244 - t253 * t160) * qJD(5) + t301) * MDP(22) - t302 * t257 + t303; (-t219 ^ 2 - t222 ^ 2) * MDP(14) + (t168 * t222 + t169 * t219 + t284) * MDP(15) + (t153 + t289) * MDP(21) + (t152 - t288) * MDP(22); (t263 * t296 + t294) * MDP(21) + ((-t161 * t244 - t155) * t251 + (-t160 * t296 - t156) * t253 + t301) * MDP(22) + t303;];
tauc = t1;
