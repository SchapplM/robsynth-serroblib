% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:15
% EndTime: 2019-12-31 18:49:20
% DurationCPUTime: 2.00s
% Computational Cost: add. (2304->230), mult. (6148->293), div. (0->0), fcn. (4635->6), ass. (0->107)
t240 = cos(pkin(8));
t243 = cos(qJ(3));
t270 = qJD(1) * qJD(3);
t265 = t243 * t270;
t229 = t240 * t265;
t239 = sin(pkin(8));
t242 = sin(qJ(3));
t266 = t242 * t270;
t213 = -t239 * t266 + t229;
t283 = t240 * t243;
t284 = t239 * t242;
t258 = -t283 + t284;
t217 = t258 * qJD(1);
t221 = t239 * t243 + t240 * t242;
t218 = t221 * qJD(1);
t241 = sin(qJ(4));
t297 = cos(qJ(4));
t267 = qJD(4) * t297;
t274 = qJD(4) * t241;
t278 = t239 * t265 + t240 * t266;
t172 = -t297 * t213 + t217 * t267 + t218 * t274 + t241 * t278;
t195 = t297 * t217 + t218 * t241;
t238 = qJD(3) + qJD(4);
t287 = t195 * t238;
t160 = -t172 + t287;
t255 = -t241 * t217 + t297 * t218;
t173 = t255 * qJD(4) + t241 * t213 + t297 * t278;
t288 = t195 ^ 2;
t289 = t255 * t238;
t304 = t255 ^ 2;
t309 = t195 * t255;
t310 = t160 * MDP(17) + MDP(15) * t309 + (-t173 + t289) * MDP(18) + (-t288 + t304) * MDP(16);
t233 = -t240 * pkin(2) - pkin(1);
t224 = t233 * qJD(1) + qJD(2);
t202 = pkin(3) * t217 + t224;
t165 = pkin(4) * t195 - qJ(5) * t255 + t202;
t307 = t165 * t195;
t306 = t195 * t202;
t174 = pkin(4) * t255 + qJ(5) * t195;
t295 = pkin(6) + qJ(2);
t225 = t295 * t239;
t222 = qJD(1) * t225;
t226 = t295 * t240;
t223 = qJD(1) * t226;
t302 = -t243 * t222 - t223 * t242;
t301 = -pkin(7) * t221 - t242 * t226;
t300 = (t239 ^ 2 + t240 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t181 = -t278 * pkin(7) - qJD(2) * t217 + t302 * qJD(3);
t252 = t221 * qJD(2);
t249 = qJD(1) * t252;
t260 = t222 * t242 - t223 * t243;
t182 = -pkin(7) * t213 + t260 * qJD(3) - t249;
t190 = -pkin(7) * t218 + t302;
t189 = qJD(3) * pkin(3) + t190;
t191 = -pkin(7) * t217 - t260;
t155 = t241 * t181 - t297 * t182 + t189 * t274 + t191 * t267;
t253 = t165 * t255 + t155;
t299 = -t255 * t202 - t155;
t276 = qJD(2) * t239;
t285 = t225 * t243;
t279 = qJD(2) * t283 - qJD(3) * t285;
t183 = t301 * qJD(3) - t242 * t276 + t279;
t219 = t258 * qJD(3);
t259 = t242 * t225 - t243 * t226;
t245 = t259 * qJD(3) - t252;
t184 = pkin(7) * t219 + t245;
t192 = -t285 + t301;
t193 = -t258 * pkin(7) - t259;
t256 = t297 * t192 - t241 * t193;
t157 = t256 * qJD(4) + t297 * t183 + t241 * t184;
t293 = t157 * t238;
t171 = t241 * t192 + t297 * t193;
t158 = t171 * qJD(4) + t241 * t183 - t297 * t184;
t292 = t158 * t238;
t269 = t297 * t191;
t164 = t241 * t189 + t269;
t291 = t164 * t238;
t282 = t241 * t191;
t167 = t297 * t190 - t282;
t280 = -pkin(3) * t267 - qJD(5) + t167;
t275 = qJD(3) * t218;
t163 = t297 * t189 - t282;
t272 = qJD(5) - t163;
t271 = qJD(1) * qJD(2);
t268 = qJD(1) * t284;
t262 = t297 * t181 + t241 * t182 + t189 * t267 - t191 * t274;
t166 = t241 * t190 + t269;
t261 = pkin(3) * t274 - t166;
t235 = t238 * qJD(5);
t154 = t235 + t262;
t254 = t163 * t238 - t262;
t251 = t221 * qJD(3);
t250 = t297 * t258;
t206 = t258 * pkin(3) + t233;
t201 = t297 * t221 - t241 * t258;
t156 = t278 * pkin(3) + t173 * pkin(4) + t172 * qJ(5) - qJD(5) * t255;
t234 = -t297 * pkin(3) - pkin(4);
t232 = pkin(3) * t241 + qJ(5);
t200 = t221 * t241 + t250;
t180 = t201 * qJD(4) - t241 * t219 + t297 * t251;
t179 = qJD(4) * t250 + t297 * t219 + t221 * t274 + t241 * t251;
t169 = t200 * pkin(4) - t201 * qJ(5) + t206;
t168 = pkin(3) * t218 + t174;
t162 = t238 * qJ(5) + t164;
t161 = -t238 * pkin(4) + t272;
t159 = pkin(3) * t251 + t180 * pkin(4) + t179 * qJ(5) - t201 * qJD(5);
t1 = [(t213 * t221 - t218 * t219) * MDP(8) + (-t213 * t258 + t219 * t217 - t218 * t251 - t221 * t278) * MDP(9) - t219 * qJD(3) * MDP(10) - t221 * qJD(3) ^ 2 * MDP(11) + (t233 * t278 + (t224 * t221 + t245) * qJD(3)) * MDP(13) + (t233 * t213 - t224 * t219 - ((-qJD(3) * t226 - t276) * t242 + t279) * qJD(3)) * MDP(14) + (-t172 * t201 - t179 * t255) * MDP(15) + (t172 * t200 - t173 * t201 + t179 * t195 - t180 * t255) * MDP(16) + (-t292 + t206 * t173 + t202 * t180 + (t195 * t251 + t278 * t200) * pkin(3)) * MDP(20) + (-t293 - t206 * t172 - t202 * t179 + (t278 * t201 + t251 * t255) * pkin(3)) * MDP(21) + (t156 * t200 + t159 * t195 + t165 * t180 + t169 * t173 - t292) * MDP(22) + (-t154 * t200 + t155 * t201 - t157 * t195 + t158 * t255 - t161 * t179 - t162 * t180 - t171 * t173 + t172 * t256) * MDP(23) + (-t156 * t201 - t159 * t255 + t165 * t179 + t169 * t172 + t293) * MDP(24) + (t154 * t171 - t155 * t256 + t156 * t169 + t157 * t162 + t158 * t161 + t159 * t165) * MDP(25) + 0.2e1 * t271 * t300 + (-t179 * MDP(17) - t180 * MDP(18)) * t238; (t275 + t278) * MDP(13) + (t229 + (-t217 - t268) * qJD(3)) * MDP(14) + (-t288 - t304) * MDP(23) + (-t161 * t255 + t162 * t195 + t156) * MDP(25) - qJD(1) ^ 2 * t300 + (-MDP(21) + MDP(24)) * (t172 + t287) + (MDP(20) + MDP(22)) * (t173 + t289); t218 * t217 * MDP(8) + (-t217 ^ 2 + t218 ^ 2) * MDP(9) + (t229 + (t217 - t268) * qJD(3)) * MDP(10) + (t275 - t278) * MDP(11) + (-t224 * t218 - t249) * MDP(13) + (t224 * t217 + t258 * t271) * MDP(14) + (t166 * t238 + (-t195 * t218 - t238 * t274) * pkin(3) + t299) * MDP(20) + (t167 * t238 + t306 + (-t218 * t255 - t238 * t267) * pkin(3) - t262) * MDP(21) + (-t168 * t195 - t261 * t238 - t253) * MDP(22) + (-t172 * t234 - t173 * t232 + (t162 + t261) * t255 + (t161 + t280) * t195) * MDP(23) + (t168 * t255 - t280 * t238 + t154 - t307) * MDP(24) + (t154 * t232 + t155 * t234 + t261 * t161 - t280 * t162 - t165 * t168) * MDP(25) + t310; (t291 + t299) * MDP(20) + (t254 + t306) * MDP(21) + (-t174 * t195 - t253 + t291) * MDP(22) + (pkin(4) * t172 - qJ(5) * t173 + (t162 - t164) * t255 + (t161 - t272) * t195) * MDP(23) + (t174 * t255 + 0.2e1 * t235 - t254 - t307) * MDP(24) + (-pkin(4) * t155 + qJ(5) * t154 - t161 * t164 + t272 * t162 - t165 * t174) * MDP(25) + t310; MDP(22) * t309 + t160 * MDP(23) + (-t238 ^ 2 - t304) * MDP(24) + (-t162 * t238 + t253) * MDP(25);];
tauc = t1;
