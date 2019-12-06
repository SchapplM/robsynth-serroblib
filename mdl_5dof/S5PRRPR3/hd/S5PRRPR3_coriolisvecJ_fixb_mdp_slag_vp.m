% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:44
% EndTime: 2019-12-05 16:19:48
% DurationCPUTime: 1.18s
% Computational Cost: add. (901->167), mult. (2327->246), div. (0->0), fcn. (1653->6), ass. (0->94)
t221 = sin(pkin(9));
t222 = cos(pkin(9));
t224 = sin(qJ(3));
t226 = cos(qJ(3));
t206 = -t221 * t224 + t222 * t226;
t201 = t206 * qJD(2);
t225 = cos(qJ(5));
t192 = t225 * t201;
t207 = t221 * t226 + t222 * t224;
t202 = t207 * qJD(3);
t195 = qJD(2) * t202;
t243 = qJD(2) * qJD(3);
t238 = t226 * t243;
t239 = t224 * t243;
t196 = -t221 * t239 + t222 * t238;
t203 = t207 * qJD(2);
t223 = sin(qJ(5));
t247 = qJD(5) * t223;
t135 = qJD(5) * t192 - t223 * t195 + t225 * t196 - t203 * t247;
t231 = t201 * t223 + t225 * t203;
t136 = qJD(5) * t231 + t225 * t195 + t196 * t223;
t162 = -t203 * t223 + t192;
t218 = qJD(3) + qJD(5);
t253 = t162 * t218;
t254 = t231 * t218;
t270 = (-t136 + t254) * MDP(17) + (-t162 ^ 2 + t231 ^ 2) * MDP(15) - t162 * MDP(14) * t231 + (t135 - t253) * MDP(16);
t256 = -qJ(4) - pkin(6);
t212 = t256 * t224;
t197 = t226 * qJD(1) + qJD(2) * t212;
t255 = qJD(3) * pkin(3);
t191 = t197 + t255;
t213 = t256 * t226;
t199 = qJD(1) * t224 - qJD(2) * t213;
t252 = t222 * t199;
t152 = t221 * t191 + t252;
t258 = pkin(7) * t201;
t144 = t152 + t258;
t240 = -pkin(3) * t226 - pkin(2);
t233 = t240 * qJD(2);
t211 = qJD(4) + t233;
t169 = -pkin(4) * t201 + t211;
t269 = t144 * t247 - t169 * t162;
t266 = -0.2e1 * t243;
t265 = MDP(5) * t224;
t264 = MDP(6) * (t224 ^ 2 - t226 ^ 2);
t263 = t224 * MDP(10) + t226 * MDP(11);
t262 = qJD(5) - t218;
t237 = qJD(3) * t256;
t198 = qJD(4) * t226 + t224 * t237;
t242 = qJD(3) * qJD(1);
t171 = qJD(2) * t198 + t226 * t242;
t200 = -qJD(4) * t224 + t226 * t237;
t172 = qJD(2) * t200 - t224 * t242;
t145 = -t171 * t221 + t222 * t172;
t138 = -pkin(7) * t196 + t145;
t146 = t222 * t171 + t221 * t172;
t139 = -pkin(7) * t195 + t146;
t261 = t225 * t138 - t223 * t139 - t169 * t231;
t259 = pkin(3) * t221;
t257 = pkin(7) * t203;
t184 = t221 * t199;
t227 = qJD(3) ^ 2;
t251 = t224 * t227;
t250 = t226 * t227;
t155 = t222 * t197 - t184;
t156 = t222 * t198 + t221 * t200;
t174 = t221 * t212 - t222 * t213;
t245 = t224 * qJD(2);
t241 = t224 * t255;
t236 = pkin(2) * t266;
t151 = t222 * t191 - t184;
t153 = -t197 * t221 - t252;
t154 = -t198 * t221 + t222 * t200;
t173 = t222 * t212 + t213 * t221;
t143 = qJD(3) * pkin(4) + t151 - t257;
t232 = -t223 * t143 - t225 * t144;
t230 = t206 * t225 - t207 * t223;
t165 = t206 * t223 + t207 * t225;
t216 = pkin(3) * t222 + pkin(4);
t215 = pkin(3) * t239;
t205 = t206 * qJD(3);
t183 = -pkin(4) * t206 + t240;
t176 = pkin(4) * t202 + t241;
t175 = pkin(3) * t245 + pkin(4) * t203;
t170 = pkin(4) * t195 + t215;
t158 = pkin(7) * t206 + t174;
t157 = -pkin(7) * t207 + t173;
t150 = -pkin(7) * t202 + t156;
t149 = t155 - t257;
t148 = -pkin(7) * t205 + t154;
t147 = t153 - t258;
t141 = qJD(5) * t165 + t225 * t202 + t205 * t223;
t140 = qJD(5) * t230 - t202 * t223 + t205 * t225;
t1 = [(-t195 * t207 - t196 * t206 + t201 * t205 + t202 * t203) * MDP(12) + (t145 * t206 + t146 * t207 - t151 * t202 + t152 * t205) * MDP(13) - t263 * t227 + (-t141 * MDP(19) - t140 * MDP(20)) * t218; 0.2e1 * t238 * t265 + t264 * t266 + MDP(7) * t250 - MDP(8) * t251 + (-pkin(6) * t250 + t224 * t236) * MDP(10) + (pkin(6) * t251 + t226 * t236) * MDP(11) + (-t145 * t207 + t146 * t206 - t151 * t205 - t152 * t202 - t154 * t203 + t156 * t201 - t173 * t196 - t174 * t195) * MDP(12) + (t145 * t173 + t146 * t174 + t151 * t154 + t152 * t156 + (t211 + t233) * t241) * MDP(13) + (t135 * t165 + t140 * t231) * MDP(14) + (t135 * t230 - t136 * t165 + t140 * t162 - t141 * t231) * MDP(15) + (t183 * t136 + t169 * t141 - t162 * t176 - t170 * t230) * MDP(19) + (t183 * t135 + t169 * t140 + t170 * t165 + t176 * t231) * MDP(20) + (t140 * MDP(16) - t141 * MDP(17) + (t148 * t225 - t150 * t223 + (-t157 * t223 - t158 * t225) * qJD(5)) * MDP(19) + (-t148 * t223 - t150 * t225 - (t157 * t225 - t158 * t223) * qJD(5)) * MDP(20)) * t218; ((t152 + t153) * t203 + (t151 - t155) * t201 + (-t195 * t221 - t196 * t222) * pkin(3)) * MDP(12) + (-t151 * t153 - t152 * t155 + (t145 * t222 + t146 * t221 - t211 * t245) * pkin(3)) * MDP(13) + (t175 * t162 - (t147 * t225 - t149 * t223) * t218 + ((-t216 * t223 - t225 * t259) * t218 + t232) * qJD(5) + t261) * MDP(19) + (-t225 * t139 - t223 * t138 - t175 * t231 + (t147 * t223 + t149 * t225) * t218 + (-(t216 * t225 - t223 * t259) * t218 - t225 * t143) * qJD(5) + t269) * MDP(20) + (t263 * pkin(2) - t226 * t265 + t264) * qJD(2) ^ 2 + t270; (-t201 ^ 2 - t203 ^ 2) * MDP(12) + (t151 * t203 - t152 * t201 + t215) * MDP(13) + (t136 + t254) * MDP(19) + (t135 + t253) * MDP(20); (t262 * t232 + t261) * MDP(19) + ((-t144 * t218 - t138) * t223 + (-t262 * t143 - t139) * t225 + t269) * MDP(20) + t270;];
tauc = t1;
