% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:20
% EndTime: 2019-12-31 18:04:24
% DurationCPUTime: 1.27s
% Computational Cost: add. (764->170), mult. (2004->248), div. (0->0), fcn. (1484->6), ass. (0->90)
t208 = sin(pkin(8));
t209 = cos(pkin(8));
t212 = sin(qJ(4));
t214 = cos(qJ(4));
t183 = t208 * t212 + t209 * t214;
t174 = t183 * qJD(4);
t170 = qJD(1) * t174;
t250 = qJD(1) * t209;
t236 = t212 * t250;
t190 = qJD(4) * t236;
t246 = qJD(4) * t214;
t237 = t208 * t246;
t171 = qJD(1) * t237 - t190;
t211 = sin(qJ(5));
t213 = cos(qJ(5));
t251 = qJD(1) * t208;
t238 = t214 * t251;
t178 = -t236 + t238;
t220 = qJD(1) * t183;
t227 = t213 * t178 - t211 * t220;
t137 = t227 * qJD(5) - t211 * t170 + t213 * t171;
t253 = t213 * t220;
t153 = t178 * t211 + t253;
t245 = qJD(5) * t211;
t221 = qJD(5) * t253 + t213 * t170 + t211 * t171 + t178 * t245;
t207 = qJD(4) + qJD(5);
t267 = t227 * t207;
t268 = t153 * t207;
t273 = (-t137 + t267) * MDP(22) + t153 * t227 * MDP(19) + (-t153 ^ 2 + t227 ^ 2) * MDP(20) + (-t221 + t268) * MDP(21);
t248 = qJD(2) * t214;
t249 = qJD(2) * t212;
t272 = t208 * t248 - t209 * t249;
t224 = t272 * qJD(1);
t191 = qJ(2) * t251 + qJD(3);
t181 = -pkin(6) * t251 + t191;
t258 = -pkin(6) + qJ(2);
t188 = t258 * t209;
t185 = qJD(1) * t188;
t226 = -t181 * t212 - t185 * t214;
t140 = pkin(7) * t170 + t226 * qJD(4) + t224;
t148 = -pkin(7) * t220 - t226;
t173 = -qJD(1) * pkin(1) - pkin(2) * t250 - qJ(3) * t251 + qJD(2);
t162 = pkin(3) * t250 - t173;
t151 = pkin(4) * t220 + t162;
t271 = t151 * t153 + t148 * t245 + (-t148 * t207 - t140) * t211;
t266 = MDP(6) + MDP(9);
t265 = t214 * t181 - t185 * t212;
t264 = 0.2e1 * t220;
t263 = qJD(5) - t207;
t139 = -pkin(7) * t171 + qJD(2) * t220 + t265 * qJD(4);
t262 = -t211 * t139 + t213 * t140 - t151 * t227;
t259 = pkin(4) * t178;
t255 = t209 * MDP(8);
t254 = t213 * t148;
t205 = t208 ^ 2;
t206 = t209 ^ 2;
t252 = t205 + t206;
t247 = qJD(4) * t212;
t243 = t205 * MDP(10);
t242 = t208 * qJD(3);
t241 = qJD(1) * qJD(2);
t240 = qJD(1) * qJD(3);
t239 = -t209 * pkin(2) - t208 * qJ(3) - pkin(1);
t234 = qJ(2) * t241;
t233 = t208 * t240;
t147 = -pkin(7) * t178 + t265;
t146 = qJD(4) * pkin(4) + t147;
t232 = -pkin(4) * t207 - t146;
t179 = t209 * pkin(3) - t239;
t229 = t206 * t234;
t184 = t208 * t214 - t209 * t212;
t157 = t183 * t213 + t184 * t211;
t158 = -t183 * t211 + t184 * t213;
t187 = t258 * t208;
t225 = -t187 * t212 - t188 * t214;
t222 = t187 * t246 - t188 * t247 + t208 * t249 + t209 * t248;
t218 = t225 * qJD(4) + t272;
t216 = qJD(1) ^ 2;
t192 = t205 * t234;
t175 = -t209 * t247 + t237;
t161 = pkin(4) * t175 + t242;
t160 = pkin(4) * t171 + t233;
t159 = pkin(4) * t183 + t179;
t150 = -pkin(7) * t183 - t225;
t149 = -pkin(7) * t184 + t187 * t214 - t188 * t212;
t144 = pkin(7) * t174 + t218;
t143 = -pkin(7) * t175 + t222;
t142 = t158 * qJD(5) - t174 * t211 + t213 * t175;
t141 = -t157 * qJD(5) - t174 * t213 - t175 * t211;
t1 = [0.2e1 * (t192 + t229) * MDP(7) + 0.2e1 * t233 * t255 + 0.2e1 * t240 * t243 + (0.2e1 * t229 + t192 + (t191 * qJD(2) + (-qJD(1) * t239 - t173) * qJD(3)) * t208) * MDP(11) + (-t170 * t184 - t174 * t178) * MDP(12) + (t170 * t183 - t171 * t184 + t174 * t220 - t175 * t178) * MDP(13) + (t162 * t175 + t179 * t171 + t264 * t242) * MDP(17) + (-t179 * t170 - t162 * t174 + (qJD(1) * t184 + t178) * t242) * MDP(18) + (t141 * t227 - t158 * t221) * MDP(19) + (-t137 * t158 - t141 * t153 - t142 * t227 + t157 * t221) * MDP(20) + (t159 * t137 + t151 * t142 + t161 * t153 + t160 * t157) * MDP(24) + (t151 * t141 + t160 * t158 - t159 * t221 + t161 * t227) * MDP(25) + 0.2e1 * t266 * t252 * t241 + (t141 * MDP(21) - t142 * MDP(22) + (-t143 * t211 + t213 * t144 + (-t149 * t211 - t150 * t213) * qJD(5)) * MDP(24) + (-t143 * t213 - t211 * t144 - (t149 * t213 - t150 * t211) * qJD(5)) * MDP(25)) * t207 + (-t174 * MDP(14) - t175 * MDP(15) + t218 * MDP(17) - t222 * MDP(18)) * qJD(4); (-qJ(2) * t206 * t216 + (-qJD(3) - t191) * t251) * MDP(11) + (t190 + (-t178 - t238) * qJD(4)) * MDP(17) + t264 * qJD(4) * MDP(18) + (-t137 - t267) * MDP(24) + (t221 + t268) * MDP(25) - (qJ(2) * MDP(7) + t266) * t252 * t216; -t216 * t243 + (-MDP(17) * t212 - MDP(18) * t214) * qJD(4) ^ 2 + (-t216 * t255 + ((qJD(2) + t173) * MDP(11) - t220 * MDP(17) - t178 * MDP(18) - t153 * MDP(24) - t227 * MDP(25)) * qJD(1)) * t208 + t207 ^ 2 * ((-t211 * t214 - t212 * t213) * MDP(24) - (-t211 * t212 + t213 * t214) * MDP(25)); t178 ^ 2 * MDP(13) + (t190 + (t178 - t238) * qJD(4)) * MDP(15) + (-t162 * t178 + t224) * MDP(17) + (-t153 * t259 - (-t147 * t211 - t254) * t207 + (t232 * t211 - t254) * qJD(5) + t262) * MDP(24) + (-t227 * t259 + (t232 * qJD(5) + t147 * t207 - t139) * t213 + t271) * MDP(25) - (-t178 * MDP(12) + (qJD(2) - t162) * MDP(18) + MDP(13) * t220) * t220 + t273; (t263 * (-t146 * t211 - t254) + t262) * MDP(24) + ((-t263 * t146 - t139) * t213 + t271) * MDP(25) + t273;];
tauc = t1;
