% Calculate Coriolis joint torque vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:37
% EndTime: 2021-01-15 15:32:41
% DurationCPUTime: 1.19s
% Computational Cost: add. (1066->195), mult. (2703->267), div. (0->0), fcn. (1792->6), ass. (0->92)
t215 = sin(qJ(3));
t217 = cos(qJ(3));
t269 = qJ(4) + pkin(6);
t235 = qJD(3) * t269;
t182 = qJD(4) * t217 - t215 * t235;
t213 = sin(pkin(8));
t214 = cos(pkin(8));
t224 = -qJD(4) * t215 - t217 * t235;
t193 = t213 * t215 - t214 * t217;
t218 = cos(qJ(2));
t227 = t193 * t218;
t261 = qJD(1) * t227 + t214 * t182 + t213 * t224;
t194 = t213 * t217 + t214 * t215;
t188 = t194 * qJD(2);
t275 = MDP(14) + MDP(17);
t252 = qJD(2) * t217;
t238 = t214 * t252;
t254 = qJD(2) * t215;
t185 = t213 * t254 - t238;
t210 = -pkin(3) * t217 - pkin(2);
t255 = qJD(1) * t218;
t191 = t210 * qJD(2) + qJD(4) - t255;
t248 = t191 * MDP(15);
t274 = (t185 * MDP(12) + t188 * MDP(13) + t248) * pkin(3);
t273 = (t215 ^ 2 - t217 ^ 2) * MDP(6);
t262 = t182 * t213 - t194 * t255 - t214 * t224;
t249 = qJD(3) * t217;
t250 = qJD(3) * t215;
t271 = t213 * t250 - t214 * t249;
t270 = -MDP(10) * t217 + MDP(11) * t215;
t187 = t194 * qJD(3);
t243 = MDP(12) + MDP(16);
t242 = MDP(13) - MDP(18);
t184 = t188 ^ 2;
t268 = qJD(2) * pkin(2);
t216 = sin(qJ(2));
t256 = qJD(1) * t216;
t201 = qJD(2) * pkin(6) + t256;
t231 = qJD(4) + t255;
t164 = -t201 * t250 + (-qJ(4) * t250 + t231 * t217) * qJD(2);
t221 = -t201 * t249 + (-qJ(4) * t249 - t231 * t215) * qJD(2);
t140 = t164 * t213 - t214 * t221;
t199 = t269 * t217;
t237 = t269 * t215;
t166 = t199 * t213 + t214 * t237;
t267 = t140 * t166;
t178 = t194 * t216;
t266 = t140 * t178;
t234 = qJ(4) * qJD(2) + t201;
t181 = t234 * t217;
t265 = t181 * t213;
t169 = t214 * t181;
t141 = t214 * t164 + t213 * t221;
t180 = t234 * t215;
t171 = qJD(3) * pkin(3) - t180;
t150 = t213 * t171 + t169;
t245 = qJD(2) * qJD(3);
t236 = t215 * t245;
t253 = qJD(2) * t216;
t192 = pkin(3) * t236 + qJD(1) * t253;
t247 = t217 * MDP(11);
t153 = -t180 * t214 - t265;
t246 = qJD(5) - t153;
t244 = qJD(3) * qJD(5);
t241 = t185 * t256;
t145 = pkin(3) * t250 + pkin(4) * t187 + qJ(5) * t271 - qJD(5) * t194;
t233 = -t145 + t256;
t230 = qJD(3) * t153 - t141;
t149 = t171 * t214 - t265;
t229 = -0.2e1 * t268;
t176 = qJD(2) * t187;
t200 = t213 * t236;
t177 = qJD(3) * t238 - t200;
t226 = pkin(4) * t176 - qJ(5) * t177 + t192;
t167 = t214 * t199 - t213 * t237;
t222 = t140 * t194 + t166 * t177 - t167 * t176 - t261 * t185;
t220 = qJD(2) ^ 2;
t219 = qJD(3) ^ 2;
t208 = -pkin(3) * t214 - pkin(4);
t206 = pkin(3) * t213 + qJ(5);
t179 = t193 * t216;
t165 = pkin(4) * t193 - qJ(5) * t194 + t210;
t156 = pkin(3) * t254 + pkin(4) * t188 + qJ(5) * t185;
t155 = -qJD(2) * t227 - t216 * t187;
t154 = -t218 * t188 + t271 * t216;
t152 = -t180 * t213 + t169;
t148 = qJD(3) * qJ(5) + t150;
t147 = pkin(4) * t185 - qJ(5) * t188 + t191;
t146 = -qJD(3) * pkin(4) + qJD(5) - t149;
t142 = -qJD(5) * t188 + t226;
t139 = t244 + t141;
t1 = [(-t141 * t179 + t149 * t154 + t150 * t155 + t266) * MDP(15) + (-t139 * t179 - t146 * t154 + t148 * t155 + t266) * MDP(19) + t243 * t185 * t253 + (t243 * t154 - t242 * t155) * qJD(3) + (-t192 * MDP(15) - t142 * MDP(19) - t220 * MDP(4) - t242 * t177 - t243 * t176 + 0.2e1 * (-t215 * MDP(10) - t247) * t245) * t218 + (-t220 * MDP(3) + (t147 * MDP(19) + t242 * t188 + t248) * qJD(2) + t270 * (t219 + t220)) * t216 + t275 * (-t154 * t188 - t155 * t185 + t179 * t176 + t177 * t178); (t176 * t210 + t187 * t191 + t192 * t193 - t241) * MDP(12) + (t177 * t210 - t191 * t271 + t192 * t194) * MDP(13) + (-t141 * t193 + t149 * t271 - t150 * t187 + t222) * MDP(14) + (t141 * t167 - t262 * t149 + t261 * t150 - t191 * t256 + t192 * t210 + t267) * MDP(15) + (t142 * t193 + t145 * t185 + t147 * t187 + t165 * t176 - t241) * MDP(16) + (-t139 * t193 - t146 * t271 - t148 * t187 + t222) * MDP(17) + (-t142 * t194 + t147 * t271 - t165 * t177) * MDP(18) + (t139 * t167 + t142 * t165 + t262 * t146 - t233 * t147 + t261 * t148 + t267) * MDP(19) + (-MDP(13) * t256 + t233 * MDP(18) + t275 * t262) * t188 + (t217 * MDP(7) - t215 * MDP(8) + t270 * pkin(6)) * t219 + (-0.2e1 * qJD(2) * t273 + t229 * t247 + (t229 * MDP(10) + 0.2e1 * MDP(5) * t252 + t274) * t215 - t243 * t262 - t242 * t261) * qJD(3); t230 * MDP(13) + (t149 * t152 - t150 * t153) * MDP(15) + (-t176 * t206 + t177 * t208) * MDP(17) + (-t230 + 0.2e1 * t244) * MDP(18) + (t139 * t206 + t140 * t208 - t146 * t152 - t147 * t156 + t246 * t148) * MDP(19) + (-t191 * MDP(12) + (t150 - t152) * MDP(14) - t147 * MDP(16) + (t148 - t152) * MDP(17) + t156 * MDP(18)) * t188 + ((-t176 * t213 - t177 * t214) * MDP(14) + (-t140 * t214 + t141 * t213) * MDP(15)) * pkin(3) + (-t215 * t217 * MDP(5) + t273) * t220 + (t191 * MDP(13) + (-t149 + t153) * MDP(14) - t156 * MDP(16) + (t146 - t246) * MDP(17) - t147 * MDP(18)) * t185 + (t268 * t247 + (MDP(10) * t268 - t274) * t215) * qJD(2) + t243 * (qJD(3) * t152 - t140); (t149 * t188 + t150 * t185 + t192) * MDP(15) + (t148 * t185 + (-qJD(5) - t146) * t188 + t226) * MDP(19) + t242 * (-t200 + (-t185 + t238) * qJD(3)) + 0.2e1 * t243 * t188 * qJD(3) + t275 * (-t185 ^ 2 - t184); t188 * t185 * MDP(16) + (-t200 + (t185 + t238) * qJD(3)) * MDP(17) + (-t184 - t219) * MDP(18) + (-qJD(3) * t148 + t147 * t188 + t140) * MDP(19);];
tauc = t1;
