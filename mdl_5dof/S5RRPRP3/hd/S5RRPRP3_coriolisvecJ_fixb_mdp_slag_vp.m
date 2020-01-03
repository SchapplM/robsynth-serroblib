% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:16
% EndTime: 2019-12-31 19:51:19
% DurationCPUTime: 1.12s
% Computational Cost: add. (1526->198), mult. (2620->253), div. (0->0), fcn. (1662->6), ass. (0->96)
t236 = cos(qJ(2));
t281 = pkin(1) * qJD(1);
t242 = -t236 * t281 + qJD(3);
t231 = sin(pkin(8));
t232 = cos(pkin(8));
t266 = t231 ^ 2 + t232 ^ 2;
t235 = cos(qJ(4));
t274 = t235 * t232;
t233 = sin(qJ(4));
t278 = t231 * t233;
t210 = -t274 + t278;
t211 = t231 * t235 + t232 * t233;
t230 = qJD(1) + qJD(2);
t194 = t211 * t230;
t216 = (-pkin(7) - qJ(3)) * t231;
t226 = t232 * pkin(7);
t217 = qJ(3) * t232 + t226;
t240 = t216 * t235 - t217 * t233;
t271 = t240 * qJD(4) - t210 * t242;
t183 = t216 * t233 + t217 * t235;
t270 = t183 * qJD(4) + t211 * t242;
t260 = MDP(16) + MDP(18);
t259 = -MDP(17) + MDP(20);
t283 = t194 ^ 2;
t282 = pkin(1) * t236;
t280 = pkin(1) * qJD(2);
t234 = sin(qJ(2));
t213 = qJ(3) * t230 + t234 * t281;
t252 = pkin(7) * t230 + t213;
t185 = t252 * t232;
t279 = t185 * t233;
t276 = t232 * MDP(7);
t262 = qJD(4) * t235;
t253 = t232 * t262;
t214 = t230 * t253;
t263 = qJD(4) * t233;
t254 = t231 * t263;
t188 = t230 * t254 - t214;
t206 = t211 * qJD(4);
t189 = t230 * t206;
t256 = qJD(1) * t280;
t222 = t234 * t256;
t238 = pkin(4) * t189 + qJ(5) * t188 + t222;
t155 = -qJD(5) * t194 + t238;
t225 = -pkin(3) * t232 - pkin(2);
t191 = t225 * t230 + t242;
t255 = t230 * t278;
t192 = -t230 * t274 + t255;
t161 = t192 * pkin(4) - t194 * qJ(5) + t191;
t273 = t155 * t210 + t161 * t206;
t205 = -t253 + t254;
t272 = -t155 * t211 + t161 * t205;
t269 = t191 * t206 + t210 * t222;
t268 = -t191 * t205 + t211 * t222;
t209 = qJD(3) * t230 + t236 * t256;
t267 = t266 * t209;
t220 = t236 * t280 + qJD(3);
t224 = pkin(1) * t234 + qJ(3);
t207 = (-pkin(7) - t224) * t231;
t208 = t224 * t232 + t226;
t241 = t207 * t235 - t208 * t233;
t159 = qJD(4) * t241 - t210 * t220;
t265 = qJD(4) * t159;
t177 = t207 * t233 + t208 * t235;
t160 = qJD(4) * t177 + t211 * t220;
t264 = qJD(4) * t160;
t184 = t252 * t231;
t165 = -t184 * t235 - t279;
t261 = qJD(5) - t165;
t257 = t234 * t280;
t251 = t213 * t266;
t250 = t266 * MDP(9);
t249 = t266 * t230;
t248 = t165 + t279;
t247 = MDP(10) * t266;
t246 = t230 * t257;
t245 = (t188 * t210 - t189 * t211 + t192 * t205 - t194 * t206) * MDP(12) + (-t188 * t211 - t194 * t205) * MDP(11) + (-MDP(13) * t205 - MDP(14) * t206) * qJD(4);
t239 = t184 * t262 + t209 * t210;
t149 = (qJD(5) - t279) * qJD(4) - t239;
t150 = -t184 * t263 + t185 * t262 + t209 * t211;
t162 = -qJD(4) * pkin(4) + t261;
t166 = -t184 * t233 + t185 * t235;
t163 = qJD(4) * qJ(5) + t166;
t244 = -t149 * t210 + t150 * t211 - t162 * t205 - t163 * t206;
t237 = -t161 * t194 - t150;
t167 = pkin(4) * t206 + qJ(5) * t205 - qJD(5) * t211;
t178 = pkin(4) * t210 - qJ(5) * t211 + t225;
t218 = t231 * t222;
t215 = t225 - t282;
t212 = -pkin(2) * t230 + t242;
t190 = t192 ^ 2;
t173 = t178 - t282;
t172 = pkin(4) * t194 + qJ(5) * t192;
t171 = t214 + (t192 - t255) * qJD(4);
t164 = t167 + t257;
t1 = [(-t222 - t246) * MDP(5) + (t231 * t246 + t218) * MDP(8) + (t220 * t249 + t267) * MDP(9) + ((t212 + (-pkin(2) - t282) * qJD(1)) * t257 + (t209 * t224 + t213 * t220) * t266) * MDP(10) + (t189 * t215 + t192 * t257 - t264 + t269) * MDP(16) + (-t188 * t215 + t194 * t257 - t265 + t268) * MDP(17) + (t164 * t192 + t173 * t189 - t264 + t273) * MDP(18) + (-t159 * t192 + t160 * t194 - t177 * t189 + t188 * t241 + t244) * MDP(19) + (-t164 * t194 + t173 * t188 + t265 + t272) * MDP(20) + (t149 * t177 - t150 * t241 + t155 * t173 + t159 * t163 + t160 * t162 + t161 * t164) * MDP(21) + t245 + (t236 * MDP(6) + t234 * t276) * (-qJD(1) - t230) * t280; -t222 * MDP(5) + t218 * MDP(8) + t267 * MDP(9) + (t189 * t225 + t269) * MDP(16) + (-t188 * t225 + t268) * MDP(17) + (t167 * t192 + t178 * t189 + t273) * MDP(18) + (-t183 * t189 + t188 * t240 - t192 * t271 + t194 * t270 + t244) * MDP(19) + (-t167 * t194 + t178 * t188 + t272) * MDP(20) + (t149 * t183 - t150 * t240 + t155 * t178 + t161 * t167 + t162 * t270 + t163 * t271) * MDP(21) + t209 * qJ(3) * t247 + (MDP(9) * t249 + t213 * t247) * qJD(3) + (t259 * t271 - t260 * t270) * qJD(4) + ((-qJD(2) * MDP(6) + (MDP(6) - t250) * t230 - MDP(10) * t251) * t236 + (-qJD(2) * t276 + (-pkin(2) * qJD(2) - t212) * MDP(10) - t161 * MDP(21) + (-MDP(8) * t231 + MDP(5) + t276) * t230 + t259 * t194 - t260 * t192) * t234) * t281 + t245; (-t230 * t251 + t222) * MDP(10) + (-t190 - t283) * MDP(19) + (t163 * t192 + (-qJD(5) - t162) * t194 + t238) * MDP(21) - t230 ^ 2 * t250 + t259 * (-t214 + (t192 + t255) * qJD(4)) + 0.2e1 * t260 * t194 * qJD(4); (-t190 + t283) * MDP(12) + t171 * MDP(13) + (-t191 * t194 - t150) * MDP(16) + t239 * MDP(17) + t237 * MDP(18) + (pkin(4) * t188 - qJ(5) * t189 - (-t163 + t166) * t194) * MDP(19) + (t172 * t194 - t239) * MDP(20) + (-pkin(4) * t150 + qJ(5) * t149 - t161 * t172 - t162 * t166 + t163 * t261) * MDP(21) + (t194 * MDP(11) + t191 * MDP(17) - t172 * MDP(18) + (t162 - t261) * MDP(19) - t161 * MDP(20)) * t192 + (t248 * MDP(17) + (0.2e1 * qJD(5) - t248) * MDP(20) + t260 * t166) * qJD(4); t194 * t192 * MDP(18) + t171 * MDP(19) + (-qJD(4) ^ 2 - t283) * MDP(20) + (-qJD(4) * t163 - t237) * MDP(21);];
tauc = t1;
