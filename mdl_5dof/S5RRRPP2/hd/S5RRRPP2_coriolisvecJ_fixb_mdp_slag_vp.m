% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:52:02
% EndTime: 2019-12-31 20:52:05
% DurationCPUTime: 1.27s
% Computational Cost: add. (1059->240), mult. (1758->293), div. (0->0), fcn. (777->4), ass. (0->115)
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t262 = qJD(3) * t227;
t284 = qJ(4) * t262 + t225 * qJD(4);
t220 = qJD(1) + qJD(2);
t226 = sin(qJ(2));
t276 = pkin(1) * qJD(1);
t251 = t226 * t276;
t195 = pkin(7) * t220 + t251;
t186 = t225 * t195;
t273 = qJ(5) * t220;
t171 = t225 * t273 - t186;
t255 = qJD(4) - t171;
t280 = pkin(3) + pkin(4);
t160 = -t280 * qJD(3) + t255;
t283 = MDP(7) * t225;
t223 = t225 ^ 2;
t224 = t227 ^ 2;
t282 = (t223 - t224) * MDP(8);
t187 = t227 * t195;
t172 = -t227 * t273 + t187;
t222 = qJD(3) * qJ(4);
t167 = t172 + t222;
t215 = t225 * qJ(4);
t218 = t227 * pkin(3);
t281 = t215 + t218;
t279 = pkin(2) * t220;
t230 = qJD(3) ^ 2;
t278 = pkin(7) * t230;
t277 = pkin(7) - qJ(5);
t275 = pkin(1) * qJD(2);
t274 = qJ(4) * t227;
t211 = pkin(1) * t226 + pkin(7);
t272 = t211 * t230;
t271 = -qJ(5) + t211;
t247 = qJD(1) * t275;
t240 = t226 * t247;
t233 = t284 * t220 - t240;
t263 = qJD(3) * t225;
t246 = t220 * t263;
t153 = -t280 * t246 + t233;
t244 = pkin(2) + t215;
t228 = cos(qJ(2));
t250 = t228 * t276;
t155 = t250 + qJD(5) + (t280 * t227 + t244) * t220;
t270 = t153 * t225 + t155 * t262;
t183 = t195 * t262;
t239 = t228 * t247;
t200 = t225 * t239;
t168 = t183 + t200;
t196 = -t250 - t279;
t269 = t196 * t262 + t225 * t240;
t241 = t220 * t251;
t268 = t227 * t241 + t250 * t263;
t266 = t223 + t224;
t265 = MDP(18) * t220;
t264 = MDP(19) * t220;
t261 = qJD(5) * t225;
t260 = qJD(5) * t227;
t166 = -t250 + (-t244 - t218) * t220;
t259 = t166 * MDP(17);
t243 = -qJD(3) * pkin(3) + qJD(4);
t176 = t186 + t243;
t258 = t176 * MDP(17);
t257 = t211 * MDP(17);
t256 = t230 * MDP(10);
t254 = qJD(5) + t155;
t253 = -MDP(14) - MDP(18);
t252 = pkin(2) + t281;
t249 = t226 * t275;
t248 = t228 * t275;
t181 = pkin(3) * t263 - t284;
t212 = -pkin(1) * t228 - pkin(2);
t245 = t220 * t262;
t209 = qJ(5) * t263;
t204 = t277 * t227;
t191 = t271 * t227;
t156 = pkin(3) * t246 - t233;
t242 = -t181 * t220 - t156;
t238 = -qJD(5) + t248;
t201 = t227 * t239;
t237 = -t195 * t263 + t201;
t188 = t212 - t281;
t221 = qJD(3) * qJD(4);
t163 = t221 + t237;
t236 = t188 * t220 - t248;
t235 = t212 * t220 - t248;
t174 = t181 + t249;
t234 = -t174 * t220 - t156 - t272;
t173 = -pkin(4) * t263 - t181;
t232 = t220 * t249 + t272;
t231 = -0.2e1 * t220 * qJD(3) * t282 + 0.2e1 * t245 * t283 + t230 * t227 * MDP(9) + t167 * t263 * MDP(20) + (t163 * t227 + t168 * t225 + t176 * t262) * MDP(15);
t219 = t220 ^ 2;
t217 = t227 * pkin(4);
t216 = 0.2e1 * t221;
t205 = t220 * t209;
t203 = t277 * t225;
t193 = t225 * t241;
t190 = t271 * t225;
t189 = t217 + t252;
t184 = t196 * t263;
t182 = (pkin(3) * t225 - t274) * t220;
t180 = qJD(3) * t204 - t261;
t179 = -pkin(7) * t263 + t209 - t260;
t178 = t187 + t222;
t177 = -t188 + t217;
t170 = (-t280 * t225 + t274) * t220;
t165 = t173 - t249;
t162 = qJD(3) * t191 + t238 * t225;
t161 = -t211 * t263 + t238 * t227 + t209;
t158 = t166 * t263;
t152 = (-qJ(5) * t262 - t261) * t220 + t168;
t151 = t153 * t227;
t149 = -t220 * t260 + t163 + t205;
t1 = [t184 * MDP(12) + t269 * MDP(13) + t158 * MDP(14) + (t156 * t188 + t166 * t174) * MDP(17) + (-qJD(3) * t162 + t151) * MDP(18) + (qJD(3) * t161 + t270) * MDP(19) + (t149 * t191 + t152 * t190 + t153 * t177 + t155 * t165 + t160 * t162 + t161 * t167) * MDP(21) + ((-qJD(1) - t220) * MDP(5) * t226 + (-qJD(1) * MDP(6) + (t266 * MDP(15) - MDP(6)) * t220) * t228) * t275 + (-t256 + t232 * MDP(13) + t234 * MDP(16) + (t168 * t211 + t176 * t248) * MDP(17) + t165 * t264 + (-t162 * t220 - t152) * MDP(20) + (t235 * MDP(12) + t236 * MDP(14) + (-t177 * t220 - t155) * MDP(18) + t191 * t220 * MDP(20) + (-MDP(15) - t257) * t178) * qJD(3)) * t225 + ((-t232 - t240) * MDP(12) + t234 * MDP(14) + (t163 * t211 + t178 * t248) * MDP(17) + t165 * t265 + (-t161 * t220 - t149) * MDP(20) + (t235 * MDP(13) + (-t166 - t236) * MDP(16) + t176 * t257 + t177 * t264 + (-t190 * t220 - t160) * MDP(20)) * qJD(3)) * t227 + t231; (t184 + t268) * MDP(12) + (-t193 + t269) * MDP(13) + (t158 + t268) * MDP(14) + t193 * MDP(16) + (-t156 * t252 + t166 * t181) * MDP(17) + (-qJD(3) * t180 + t151 + t268) * MDP(18) + (qJD(3) * t179 + t193 + t270) * MDP(19) + (t149 * t204 + t152 * t203 + t153 * t189 + t155 * t173 + t160 * t180 + t167 * t179) * MDP(21) + (-t256 + t242 * MDP(16) + t173 * t264 + (-t180 * t220 - t152) * MDP(20) + (t168 * MDP(17) + (MDP(13) - MDP(16)) * t230) * pkin(7) + (-t155 * MDP(18) + (-pkin(7) * MDP(17) - MDP(15)) * t178 + (-pkin(2) * MDP(12) - MDP(14) * t252 - t189 * MDP(18) + t204 * MDP(20)) * t220) * qJD(3)) * t225 + (((-qJD(2) + t220) * MDP(5) - t259 + t155 * MDP(21)) * t226 + (-qJD(2) * MDP(6) + (-t160 * MDP(21) - t258) * t225 + (MDP(6) + (-MDP(15) + MDP(20)) * t266) * t220) * t228) * t276 + ((-t240 - t278) * MDP(12) + (t242 - t278) * MDP(14) + (pkin(7) * t163 - t178 * t250) * MDP(17) + t173 * t265 + (-t179 * t220 - t149) * MDP(20) - t167 * MDP(21) * t250 + ((t250 - t279) * MDP(13) + (t220 * t252 - t166 - t250) * MDP(16) + pkin(7) * t258 + (t189 * t220 - t250) * MDP(19) + (-t203 * t220 - t160) * MDP(20)) * qJD(3)) * t227 + t231; -t201 * MDP(13) + (t201 + t216) * MDP(16) + (-t176 * t187 - pkin(3) * t168 + qJ(4) * t163 - t166 * t182 + (qJD(4) + t186) * t178) * MDP(17) + (qJD(3) * t172 - t183) * MDP(18) + (-qJD(3) * t171 + t205 + t216 + t237) * MDP(19) + (qJ(4) * t149 - t152 * t280 - t155 * t170 - t160 * t172 + t255 * t167) * MDP(21) + (-MDP(12) + t253) * t200 + (-t227 * t283 + t282) * t219 + ((-t196 * MDP(12) - t166 * MDP(14) + (t178 - t222) * MDP(15) + t182 * MDP(16) + t254 * MDP(18) - t170 * MDP(19)) * t225 + (-t196 * MDP(13) + t182 * MDP(14) + (-t176 + t243) * MDP(15) + t166 * MDP(16) + (qJ(5) * qJD(3) - t170) * MDP(18) - t254 * MDP(19)) * t227) * t220; (-qJD(3) * t178 + t168) * MDP(17) + (-qJ(5) * t245 - qJD(3) * t167 + t168) * MDP(21) + (MDP(16) + MDP(19)) * (-t219 * t223 - t230) + (t253 * t227 * t219 + (-t254 * MDP(21) + t259) * t220) * t225; t233 * MDP(21) - t266 * MDP(20) * t219 + ((t160 * t225 + t167 * t227) * MDP(21) + (0.2e1 * t227 * MDP(19) + (-t280 * MDP(21) - 0.2e1 * MDP(18)) * t225) * qJD(3)) * t220;];
tauc = t1;
