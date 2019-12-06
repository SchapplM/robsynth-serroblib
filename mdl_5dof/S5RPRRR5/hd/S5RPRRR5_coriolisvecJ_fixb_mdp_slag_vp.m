% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:31
% EndTime: 2019-12-05 18:16:35
% DurationCPUTime: 0.87s
% Computational Cost: add. (893->141), mult. (1721->206), div. (0->0), fcn. (1058->8), ass. (0->93)
t272 = pkin(7) + pkin(8);
t211 = sin(qJ(5));
t212 = sin(qJ(4));
t214 = cos(qJ(5));
t215 = cos(qJ(4));
t184 = t211 * t215 + t214 * t212;
t206 = qJD(1) + qJD(3);
t178 = t184 * t206;
t205 = qJD(4) + qJD(5);
t267 = qJD(5) - t205;
t271 = MDP(8) * t212;
t270 = (t212 ^ 2 - t215 ^ 2) * MDP(9);
t199 = cos(pkin(9)) * pkin(1) + pkin(2);
t213 = sin(qJ(3));
t216 = cos(qJ(3));
t265 = pkin(1) * sin(pkin(9));
t269 = t216 * t199 - t213 * t265;
t191 = t199 * qJD(1);
t240 = qJD(1) * t265;
t172 = t213 * t191 + t216 * t240;
t233 = t206 * t272 + t172;
t152 = t215 * qJD(2) - t233 * t212;
t268 = t212 * MDP(13) + t215 * MDP(14);
t153 = t212 * qJD(2) + t233 * t215;
t264 = t206 * pkin(3);
t263 = t215 * pkin(4);
t222 = t213 * t199 + t216 * t265;
t180 = pkin(7) + t222;
t262 = -pkin(8) - t180;
t171 = t216 * t191 - t213 * t240;
t261 = t171 * t205;
t260 = t172 * t206;
t175 = t222 * qJD(3);
t259 = t175 * t206;
t258 = t211 * t212;
t256 = t214 * t153;
t254 = t214 * t215;
t201 = -pkin(3) - t263;
t156 = t201 * t206 - t171;
t230 = qJD(3) * t240;
t248 = qJD(3) * t191;
t169 = t213 * t248 + t216 * t230;
t245 = t212 * qJD(4);
t239 = pkin(4) * t245;
t157 = t206 * t239 + t169;
t183 = -t254 + t258;
t160 = t205 * t183;
t252 = -t156 * t160 + t157 * t184;
t161 = t205 * t184;
t251 = t156 * t161 + t157 * t183;
t163 = -t171 - t264;
t243 = t215 * qJD(4);
t250 = t163 * t243 + t169 * t212;
t217 = qJD(4) ^ 2;
t242 = t217 * MDP(11);
t241 = pkin(4) * t206 * t212;
t238 = t206 * t258;
t237 = t206 * t254;
t236 = qJD(4) * t272;
t235 = t206 * t243;
t149 = qJD(4) * pkin(4) + t152;
t234 = -pkin(4) * t205 - t149;
t232 = qJD(4) * t262;
t179 = -pkin(3) - t269;
t229 = -t172 + t239;
t168 = -t213 * t230 + t216 * t248;
t227 = pkin(7) * t217 - t260;
t226 = qJD(4) * (t171 - t264);
t224 = t180 * t217 + t259;
t174 = t269 * qJD(3);
t223 = qJD(4) * (t179 * t206 - t174);
t137 = t152 * qJD(4) + t215 * t168;
t138 = -t153 * qJD(4) - t212 * t168;
t221 = -t211 * t137 + t214 * t138 - t156 * t178;
t150 = qJD(5) * t237 - t205 * t238 + t214 * t235;
t176 = -t237 + t238;
t220 = t178 * t176 * MDP(15) + (-t176 ^ 2 + t178 ^ 2) * MDP(16) + (t176 * t205 + t150) * MDP(17);
t151 = t161 * t206;
t219 = (-t150 * t183 - t184 * t151 + t160 * t176 - t178 * t161) * MDP(16) + (t150 * t184 - t178 * t160) * MDP(15) - 0.2e1 * t206 * qJD(4) * t270 + 0.2e1 * t235 * t271 + t217 * t215 * MDP(10) + (-t160 * MDP(17) - t161 * MDP(18)) * t205;
t218 = t156 * t176 + (t153 * t267 - t138) * t211;
t203 = t215 * pkin(8);
t196 = t215 * pkin(7) + t203;
t195 = t272 * t212;
t186 = t215 * t236;
t185 = t212 * t236;
t170 = t179 - t263;
t167 = t175 + t239;
t166 = t215 * t180 + t203;
t165 = t262 * t212;
t158 = t163 * t245;
t145 = -t212 * t174 + t215 * t232;
t144 = t215 * t174 + t212 * t232;
t1 = [(-t169 - t259) * MDP(6) + (-t174 * t206 - t168) * MDP(7) + t158 * MDP(13) + t250 * MDP(14) + (t167 * t176 + t170 * t151 + (-t211 * t144 + t214 * t145 + (-t165 * t211 - t166 * t214) * qJD(5)) * t205 + t251) * MDP(20) + (t167 * t178 + t170 * t150 - (t214 * t144 + t211 * t145 + (t165 * t214 - t166 * t211) * qJD(5)) * t205 + t252) * MDP(21) + ((-t169 - t224) * MDP(13) + MDP(14) * t223) * t215 + (MDP(13) * t223 + t224 * MDP(14) - t242) * t212 + t219; -t268 * t217 + (-MDP(20) * t161 + MDP(21) * t160) * t205; (-t169 + t260) * MDP(6) + (t171 * t206 - t168) * MDP(7) - t212 * t242 + (t158 + t212 * t226 + (-t169 - t227) * t215) * MDP(13) + (t227 * t212 + t215 * t226 + t250) * MDP(14) + (t201 * t151 + (t211 * t185 - t214 * t186 + (t195 * t211 - t196 * t214) * qJD(5)) * t205 + t229 * t176 + t184 * t261 + t251) * MDP(20) + (t201 * t150 - (-t214 * t185 - t211 * t186 + (-t195 * t214 - t196 * t211) * qJD(5)) * t205 + t229 * t178 - t183 * t261 + t252) * MDP(21) + t219; (-t176 * t241 - (-t211 * t152 - t256) * t205 + (t234 * t211 - t256) * qJD(5) + t221) * MDP(20) + (-t178 * t241 + (t234 * qJD(5) + t152 * t205 - t137) * t214 + t218) * MDP(21) + t220 + t268 * (-t163 * t206 - t168) + (-t215 * t271 + t270) * t206 ^ 2; (t221 + t267 * (-t211 * t149 - t256)) * MDP(20) + ((-t267 * t149 - t137) * t214 + t218) * MDP(21) + t220;];
tauc = t1;
