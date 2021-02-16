% Calculate Coriolis joint torque vector for
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:57
% EndTime: 2021-01-15 11:34:03
% DurationCPUTime: 1.52s
% Computational Cost: add. (1155->207), mult. (2561->286), div. (0->0), fcn. (1710->6), ass. (0->97)
t254 = cos(pkin(8));
t258 = cos(qJ(3));
t273 = qJD(1) * qJD(3);
t269 = t258 * t273;
t236 = t254 * t269;
t253 = sin(pkin(8));
t256 = sin(qJ(3));
t283 = qJD(1) * t256;
t270 = t253 * t283;
t214 = qJD(3) * t270 - t236;
t230 = t253 * t258 + t254 * t256;
t215 = t230 * t273;
t282 = qJD(1) * t258;
t224 = t254 * t282 - t270;
t255 = sin(qJ(5));
t257 = cos(qJ(5));
t277 = qJD(5) * t255;
t220 = t230 * qJD(1);
t289 = t257 * t220;
t153 = -qJD(5) * t289 + t255 * t214 - t257 * t215 - t224 * t277;
t264 = -t220 * t255 + t257 * t224;
t154 = t264 * qJD(5) - t257 * t214 - t215 * t255;
t178 = t224 * t255 + t289;
t248 = qJD(3) + qJD(5);
t291 = t178 * t248;
t292 = t264 * t248;
t306 = t178 * t264 * MDP(18) + (-t154 + t292) * MDP(21) + (-t178 ^ 2 + t264 ^ 2) * MDP(19) + (t153 + t291) * MDP(20);
t305 = (t258 * MDP(12) - t256 * MDP(13)) * qJ(2) - t258 * t256 * MDP(7) + (t256 ^ 2 - t258 ^ 2) * MDP(8);
t259 = -pkin(1) - pkin(6);
t237 = t259 * qJD(1) + qJD(2);
t219 = -qJ(4) * t282 + t258 * t237;
t207 = qJD(3) * pkin(3) + t219;
t218 = -qJ(4) * t283 + t237 * t256;
t290 = t254 * t218;
t170 = t253 * t207 + t290;
t295 = pkin(7) * t220;
t162 = t170 - t295;
t233 = pkin(3) * t283 + qJD(1) * qJ(2) + qJD(4);
t190 = pkin(4) * t220 + t233;
t304 = t162 * t277 + t190 * t178;
t300 = MDP(12) * t256 + MDP(13) * t258;
t299 = qJD(5) - t248;
t278 = qJD(4) * t258;
t281 = qJD(3) * t256;
t191 = -t237 * t281 + (qJ(4) * t281 - t278) * qJD(1);
t279 = qJD(4) * t256;
t280 = qJD(3) * t258;
t192 = t237 * t280 + (-qJ(4) * t280 - t279) * qJD(1);
t165 = t254 * t191 - t192 * t253;
t156 = pkin(7) * t215 + t165;
t166 = t253 * t191 + t254 * t192;
t157 = pkin(7) * t214 + t166;
t297 = t257 * t156 - t255 * t157 - t190 * t264;
t296 = pkin(3) * t253;
t294 = pkin(7) * t224;
t203 = t253 * t218;
t288 = qJ(4) - t259;
t216 = t288 * t281 - t278;
t235 = t288 * t258;
t217 = -qJD(3) * t235 - t279;
t172 = t253 * t216 + t254 * t217;
t176 = t254 * t219 - t203;
t234 = t288 * t256;
t189 = -t254 * t234 - t253 * t235;
t232 = pkin(3) * t269 + qJD(1) * qJD(2);
t242 = t256 * pkin(3) + qJ(2);
t238 = pkin(3) * t280 + qJD(2);
t272 = pkin(3) * t282;
t268 = -qJ(2) * MDP(6) - MDP(5);
t169 = t254 * t207 - t203;
t171 = t254 * t216 - t217 * t253;
t175 = -t219 * t253 - t290;
t188 = t234 * t253 - t254 * t235;
t161 = qJD(3) * pkin(4) + t169 - t294;
t265 = -t255 * t161 - t257 * t162;
t229 = -t253 * t256 + t254 * t258;
t183 = t229 * t257 - t230 * t255;
t182 = t229 * t255 + t230 * t257;
t222 = t253 * t281 - t254 * t280;
t223 = t230 * qJD(3);
t263 = t165 * t229 + t166 * t230 - t169 * t223 - t170 * t222;
t158 = -t182 * qJD(5) + t222 * t255 - t257 * t223;
t159 = t183 * qJD(5) - t257 * t222 - t223 * t255;
t261 = qJD(1) ^ 2;
t260 = qJD(3) ^ 2;
t243 = pkin(3) * t254 + pkin(4);
t208 = pkin(4) * t230 + t242;
t196 = pkin(4) * t224 + t272;
t193 = -pkin(4) * t222 + t238;
t184 = -pkin(4) * t214 + t232;
t174 = -pkin(7) * t230 + t189;
t173 = -pkin(7) * t229 + t188;
t168 = t176 - t294;
t167 = t175 + t295;
t164 = pkin(7) * t222 + t172;
t163 = pkin(7) * t223 + t171;
t1 = [(qJD(3) * t171 - t214 * t242 + t220 * t238 - t222 * t233 + t230 * t232) * MDP(14) + (-qJD(3) * t172 - t215 * t242 - t223 * t233 + t224 * t238 + t229 * t232) * MDP(15) + (-t171 * t224 - t172 * t220 + t188 * t215 + t189 * t214 - t263) * MDP(16) + (t165 * t188 + t166 * t189 + t169 * t171 + t170 * t172 + t232 * t242 + t233 * t238) * MDP(17) + (t153 * t183 + t158 * t264) * MDP(18) + (-t153 * t182 - t154 * t183 - t158 * t178 - t159 * t264) * MDP(19) + (t208 * t154 + t190 * t159 + t193 * t178 + t184 * t182) * MDP(23) + (t208 * t153 + t190 * t158 + t184 * t183 + t193 * t264) * MDP(24) + (t158 * MDP(20) - t159 * MDP(21) + (t163 * t257 - t164 * t255) * MDP(23) + (-t163 * t255 - t164 * t257) * MDP(24) + ((-t173 * t255 - t174 * t257) * MDP(23) + (-t173 * t257 + t174 * t255) * MDP(24)) * qJD(5)) * t248 + ((-MDP(13) * t259 - MDP(10)) * t258 + (-MDP(12) * t259 - MDP(9)) * t256) * t260 + (0.2e1 * (-t268 + t300) * qJD(2) + 0.2e1 * t305 * qJD(3)) * qJD(1); (-qJD(1) * t220 - qJD(3) * t223) * MDP(14) + (-qJD(1) * t224 + qJD(3) * t222) * MDP(15) + (t214 * t230 + t215 * t229 + t220 * t222 + t223 * t224) * MDP(16) + (-qJD(1) * t233 + t263) * MDP(17) + (-qJD(1) * t178 + t158 * t248) * MDP(23) + (-qJD(1) * t264 - t159 * t248) * MDP(24) + t268 * t261 + t300 * (-t260 - t261); (-qJD(3) * t175 - t220 * t272 - t224 * t233 + t165) * MDP(14) + (qJD(3) * t176 + t220 * t233 - t224 * t272 - t166) * MDP(15) + ((t170 + t175) * t224 + (-t169 + t176) * t220 + (t214 * t253 + t215 * t254) * pkin(3)) * MDP(16) + (-t169 * t175 - t170 * t176 + (t165 * t254 + t166 * t253 - t233 * t282) * pkin(3)) * MDP(17) + (-t196 * t178 - (t167 * t257 - t168 * t255) * t248 + ((-t243 * t255 - t257 * t296) * t248 + t265) * qJD(5) + t297) * MDP(23) + (-t257 * t157 - t255 * t156 - t196 * t264 + (t167 * t255 + t168 * t257) * t248 + (-(t243 * t257 - t255 * t296) * t248 - t257 * t161) * qJD(5) + t304) * MDP(24) - t305 * t261 + t306; t236 * MDP(14) + (-t220 ^ 2 - t224 ^ 2) * MDP(16) + (t169 * t224 + t170 * t220 + t232) * MDP(17) + (t154 + t292) * MDP(23) + (t153 - t291) * MDP(24) + ((t224 - t270) * MDP(14) + (-t253 * t282 - t254 * t283 - t220) * MDP(15)) * qJD(3); (t299 * t265 + t297) * MDP(23) + ((-t162 * t248 - t156) * t255 + (-t299 * t161 - t157) * t257 + t304) * MDP(24) + t306;];
tauc = t1;
