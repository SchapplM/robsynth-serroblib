% Calculate Coriolis joint torque vector for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:14
% EndTime: 2021-01-15 15:42:19
% DurationCPUTime: 1.56s
% Computational Cost: add. (987->193), mult. (2551->276), div. (0->0), fcn. (1789->6), ass. (0->93)
t232 = sin(pkin(9));
t233 = cos(pkin(9));
t237 = cos(qJ(3));
t253 = qJD(2) * t237;
t247 = t233 * t253;
t235 = sin(qJ(3));
t254 = qJD(2) * t235;
t208 = t232 * t254 - t247;
t236 = cos(qJ(5));
t199 = t236 * t208;
t215 = t232 * t237 + t233 * t235;
t210 = t215 * qJD(3);
t202 = qJD(2) * t210;
t252 = qJD(3) * t235;
t246 = qJD(2) * t252;
t222 = t232 * t246;
t203 = qJD(3) * t247 - t222;
t234 = sin(qJ(5));
t251 = qJD(5) * t234;
t256 = qJD(2) * t215;
t142 = -qJD(5) * t199 - t234 * t202 + t236 * t203 - t251 * t256;
t169 = -t234 * t256 - t199;
t241 = t208 * t234 - t236 * t256;
t240 = t241 * qJD(5) - t236 * t202 - t234 * t203;
t229 = qJD(3) + qJD(5);
t262 = t169 * t229;
t263 = t241 * t229;
t276 = t169 * MDP(16) * t241 + (-t169 ^ 2 + t241 ^ 2) * MDP(17) + (t142 - t262) * MDP(18) + (t240 - t263) * MDP(19);
t265 = -qJ(4) - pkin(6);
t220 = t265 * t235;
t204 = t237 * qJD(1) + qJD(2) * t220;
t198 = qJD(3) * pkin(3) + t204;
t221 = t265 * t237;
t206 = qJD(1) * t235 - qJD(2) * t221;
t261 = t233 * t206;
t159 = t232 * t198 + t261;
t267 = pkin(7) * t208;
t151 = t159 - t267;
t227 = -pkin(3) * t237 - pkin(2);
t255 = qJD(2) * t227;
t219 = qJD(4) + t255;
t176 = pkin(4) * t208 + t219;
t275 = t151 * t251 - t176 * t169;
t245 = qJD(3) * t265;
t205 = qJD(4) * t237 + t235 * t245;
t250 = qJD(3) * qJD(1);
t178 = t205 * qJD(2) + t237 * t250;
t207 = -qJD(4) * t235 + t237 * t245;
t179 = t207 * qJD(2) - t235 * t250;
t152 = -t178 * t232 + t233 * t179;
t145 = -pkin(7) * t203 + t152;
t153 = t233 * t178 + t232 * t179;
t146 = -pkin(7) * t202 + t153;
t274 = t236 * t145 - t234 * t146 + t176 * t241;
t272 = (t235 ^ 2 - t237 ^ 2) * MDP(6);
t258 = MDP(11) * t237;
t271 = MDP(10) * t235 + t258;
t270 = qJD(5) - t229;
t268 = pkin(3) * t232;
t266 = pkin(7) * t256;
t264 = MDP(5) * t237;
t191 = t232 * t206;
t162 = t233 * t204 - t191;
t163 = t233 * t205 + t232 * t207;
t181 = t232 * t220 - t233 * t221;
t249 = 0.2e1 * qJD(2);
t248 = pkin(3) * t254;
t158 = t233 * t198 - t191;
t160 = -t204 * t232 - t261;
t161 = -t205 * t232 + t233 * t207;
t180 = t233 * t220 + t221 * t232;
t150 = qJD(3) * pkin(4) + t158 - t266;
t242 = -t234 * t150 - t236 * t151;
t214 = t232 * t235 - t233 * t237;
t171 = t214 * t236 + t215 * t234;
t172 = -t214 * t234 + t215 * t236;
t238 = qJD(3) ^ 2;
t226 = pkin(3) * t233 + pkin(4);
t224 = pkin(3) * t246;
t213 = t214 * qJD(3);
t190 = pkin(4) * t214 + t227;
t183 = pkin(3) * t252 + pkin(4) * t210;
t182 = pkin(4) * t256 + t248;
t177 = pkin(4) * t202 + t224;
t165 = -pkin(7) * t214 + t181;
t164 = -pkin(7) * t215 + t180;
t157 = -pkin(7) * t210 + t163;
t156 = t162 - t266;
t155 = pkin(7) * t213 + t161;
t154 = t160 + t267;
t148 = t172 * qJD(5) + t236 * t210 - t213 * t234;
t147 = -t171 * qJD(5) - t210 * t234 - t213 * t236;
t1 = [(-t202 * t215 + t203 * t214 + t208 * t213 + t210 * t256) * MDP(14) + (-t152 * t214 + t153 * t215 - t158 * t210 - t159 * t213) * MDP(15) - t271 * t238 + (-t148 * MDP(21) - t147 * MDP(22)) * t229 + (-t210 * MDP(12) + t213 * MDP(13)) * qJD(3); (t202 * t227 + t210 * t219) * MDP(12) + (t203 * t227 - t213 * t219) * MDP(13) + (-t152 * t215 - t153 * t214 + t158 * t213 - t159 * t210 - t161 * t256 - t163 * t208 - t180 * t203 - t181 * t202) * MDP(14) + (t152 * t180 + t153 * t181 + t158 * t161 + t159 * t163) * MDP(15) + (t142 * t172 - t147 * t241) * MDP(16) + (-t142 * t171 + t147 * t169 + t148 * t241 + t172 * t240) * MDP(17) + (t176 * t148 - t169 * t183 + t177 * t171 - t190 * t240) * MDP(21) + (t190 * t142 + t176 * t147 + t177 * t172 - t183 * t241) * MDP(22) + (t147 * MDP(18) - t148 * MDP(19) + (t155 * t236 - t157 * t234) * MDP(21) + (-t155 * t234 - t157 * t236) * MDP(22) + ((-t164 * t234 - t165 * t236) * MDP(21) + (-t164 * t236 + t165 * t234) * MDP(22)) * qJD(5)) * t229 + (t237 * MDP(7) - t235 * MDP(8) + (-MDP(10) * t237 + MDP(11) * t235) * pkin(6)) * t238 + (t161 * MDP(12) - t163 * MDP(13) + (-pkin(2) * t258 - t272) * t249 + ((-MDP(10) * pkin(2) + t264) * t249 + ((qJD(2) * t214 + t208) * MDP(12) + 0.2e1 * t256 * MDP(13) + (t219 + t255) * MDP(15)) * pkin(3)) * t235) * qJD(3); (-qJD(3) * t160 - t208 * t248 - t219 * t256 + t152) * MDP(12) + (qJD(3) * t162 + t208 * t219 - t248 * t256 - t153) * MDP(13) + ((t159 + t160) * t256 + (-t158 + t162) * t208 + (-t202 * t232 - t203 * t233) * pkin(3)) * MDP(14) + (-t158 * t160 - t159 * t162 + (t152 * t233 + t153 * t232 - t219 * t254) * pkin(3)) * MDP(15) + (t182 * t169 - (t154 * t236 - t156 * t234) * t229 + ((-t226 * t234 - t236 * t268) * t229 + t242) * qJD(5) + t274) * MDP(21) + (-t236 * t146 - t234 * t145 + t182 * t241 + (t154 * t234 + t236 * t156) * t229 + (-(t226 * t236 - t234 * t268) * t229 - t236 * t150) * qJD(5) + t275) * MDP(22) + (pkin(2) * t271 - t235 * t264 + t272) * qJD(2) ^ 2 + t276; -t222 * MDP(13) + (-t208 ^ 2 - t256 ^ 2) * MDP(14) + (t158 * t256 + t159 * t208 + t224) * MDP(15) + (-t240 - t263) * MDP(21) + (t142 + t262) * MDP(22) + ((t232 * t253 + t233 * t254 + t256) * MDP(12) + (-t208 + t247) * MDP(13)) * qJD(3); (t242 * t270 + t274) * MDP(21) + ((-t151 * t229 - t145) * t234 + (-t150 * t270 - t146) * t236 + t275) * MDP(22) + t276;];
tauc = t1;
