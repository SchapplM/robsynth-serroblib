% Calculate Coriolis joint torque vector for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:31
% EndTime: 2021-01-15 16:23:35
% DurationCPUTime: 1.19s
% Computational Cost: add. (1127->185), mult. (2788->247), div. (0->0), fcn. (1814->4), ass. (0->103)
t201 = cos(qJ(3));
t260 = pkin(6) + pkin(7);
t184 = t260 * t201;
t200 = sin(qJ(3));
t239 = qJD(1) * t200;
t168 = qJD(2) * t184 + t239;
t199 = sin(qJ(4));
t162 = t199 * t168;
t231 = qJD(2) * t260;
t167 = t201 * qJD(1) - t200 * t231;
t256 = qJD(3) * pkin(3);
t165 = t167 + t256;
t259 = cos(qJ(4));
t222 = t259 * t165 - t162;
t179 = t199 * t201 + t200 * t259;
t238 = qJD(2) * t179;
t252 = t238 * qJ(5);
t137 = -t252 + t222;
t196 = qJD(3) + qJD(4);
t267 = t196 * MDP(19);
t234 = qJD(2) * qJD(3);
t266 = -0.2e1 * t234;
t265 = MDP(5) * t201;
t264 = MDP(6) * (t200 ^ 2 - t201 ^ 2);
t227 = t259 * qJD(4);
t263 = t259 * qJD(3) + t227;
t262 = MDP(10) * t200 + MDP(11) * t201;
t261 = t238 ^ 2;
t257 = pkin(3) * t196;
t249 = t199 * t200;
t217 = t196 * t249;
t229 = t259 * t201;
t218 = qJD(2) * t229;
t243 = t196 * t218;
t150 = qJD(2) * t217 - t243;
t255 = t150 * qJ(5);
t237 = qJD(2) * t200;
t170 = t199 * t237 - t218;
t254 = t170 * qJ(5);
t253 = t170 * t196;
t194 = -pkin(3) * t201 - pkin(2);
t182 = t194 * qJD(2);
t250 = t182 * t238;
t202 = qJD(3) ^ 2;
t248 = t200 * t202;
t247 = t201 * t202;
t136 = pkin(4) * t196 + t137;
t246 = t136 - t137;
t153 = t196 * t179;
t151 = t153 * qJD(2);
t152 = -t263 * t201 + t217;
t245 = -t179 * t151 + t152 * t170;
t244 = t259 * t167 - t162;
t236 = qJD(4) * t199;
t225 = pkin(4) * t170 + qJD(5);
t154 = t182 + t225;
t235 = qJD(5) + t154;
t233 = pkin(3) * t237;
t232 = t200 * t256;
t230 = qJD(3) * t260;
t164 = t259 * t168;
t226 = t200 * t234;
t191 = pkin(3) * t226;
t146 = pkin(4) * t151 + t191;
t224 = pkin(2) * t266;
t159 = t167 * qJD(3);
t160 = (-t201 * t231 - t239) * qJD(3);
t223 = -t199 * t159 + t259 * t160;
t221 = -t167 * t199 - t164;
t220 = t196 * t200;
t178 = -t229 + t249;
t216 = -t150 * t178 + t153 * t238;
t215 = -t199 * t165 - t164;
t183 = t260 * t200;
t214 = t199 * t183 - t184 * t259;
t169 = t170 ^ 2;
t213 = t238 * t170 * MDP(12) + (-qJD(2) * t199 * t220 + t243 + t253) * MDP(14) + (-t169 + t261) * MDP(13);
t180 = t200 * t230;
t181 = t201 * t230;
t212 = -t259 * t180 - t199 * t181 - t183 * t227 - t184 * t236;
t211 = qJD(4) * t215 + t223;
t210 = qJD(4) * t214 + t199 * t180 - t259 * t181;
t209 = t159 * t259 + t199 * t160 + t165 * t227 - t168 * t236;
t208 = t211 + t255;
t207 = t182 * t170 - t209;
t206 = -t151 * qJ(5) + t209;
t205 = (-t164 + (-t165 - t257) * t199) * qJD(4) + t223;
t204 = t170 * t235 - t206;
t193 = pkin(3) * t259 + pkin(4);
t185 = t227 * t257;
t161 = pkin(4) * t178 + t194;
t156 = pkin(4) * t238 + t233;
t149 = pkin(4) * t153 + t232;
t148 = -t178 * qJ(5) - t214;
t147 = -t179 * qJ(5) - t183 * t259 - t199 * t184;
t140 = -t252 + t244;
t139 = t221 + t254;
t138 = -t215 - t254;
t135 = t152 * qJ(5) - t179 * qJD(5) + t210;
t134 = -qJ(5) * t153 - qJD(5) * t178 + t212;
t133 = -qJD(5) * t238 + t208;
t132 = -t170 * qJD(5) + t206;
t1 = [(t216 + t245) * MDP(21) + (t132 * t179 - t133 * t178 - t136 * t153 - t138 * t152) * MDP(22) - t262 * t202 + ((-MDP(17) - MDP(19)) * t153 + (MDP(18) + MDP(20)) * t152) * t196; 0.2e1 * t226 * t265 + t264 * t266 + MDP(7) * t247 - MDP(8) * t248 + (-pkin(6) * t247 + t200 * t224) * MDP(10) + (pkin(6) * t248 + t201 * t224) * MDP(11) + (-t150 * t179 - t152 * t238) * MDP(12) + (-t216 + t245) * MDP(13) + (t194 * t151 + t182 * t153 + t170 * t232 + t178 * t191) * MDP(17) + (-t194 * t150 - t182 * t152 + 0.2e1 * t238 * t232) * MDP(18) + (t146 * t178 + t149 * t170 + t151 * t161 + t153 * t154) * MDP(19) + (t146 * t179 + t149 * t238 - t150 * t161 - t152 * t154) * MDP(20) + (-t132 * t178 - t133 * t179 - t134 * t170 - t135 * t238 + t136 * t152 - t138 * t153 + t147 * t150 - t148 * t151) * MDP(21) + (t132 * t148 + t133 * t147 + t134 * t138 + t135 * t136 + t146 * t161 + t149 * t154) * MDP(22) + (-t152 * MDP(14) - t153 * MDP(15) + MDP(17) * t210 - MDP(18) * t212 + t135 * MDP(19) - t134 * MDP(20)) * t196; (-t170 * t233 - t196 * t221 + t205 - t250) * MDP(17) + (t196 * t244 - t233 * t238 - t185 + t207) * MDP(18) + (-t139 * t196 - t156 * t170 - t235 * t238 + t205 + t255) * MDP(19) + (t140 * t196 - t156 * t238 - t185 + t204) * MDP(20) + (t193 * t150 + (t138 + t139) * t238 + (-t136 + t140) * t170 + (-t151 * t199 + (-t170 * t259 + t199 * t238) * qJD(4)) * pkin(3)) * MDP(21) + (t133 * t193 - t136 * t139 - t138 * t140 - t154 * t156 + (t132 * t199 + (-t136 * t199 + t138 * t259) * qJD(4)) * pkin(3)) * MDP(22) + t213 + (t262 * pkin(2) - t200 * t265 + t264) * qJD(2) ^ 2; (-t196 * t215 + t211 - t250) * MDP(17) + (t196 * t222 + t207) * MDP(18) + (t138 * t196 + (-t154 - t225) * t238 + t208) * MDP(19) + (-pkin(4) * t261 + t137 * t196 + t204) * MDP(20) + (pkin(4) * t150 - t170 * t246) * MDP(21) + (t246 * t138 + (-t154 * t238 + t133) * pkin(4)) * MDP(22) + t213; t238 * t267 + (t243 - t253) * MDP(20) + (-t169 - t261) * MDP(21) + (t136 * t238 + t138 * t170 + t146) * MDP(22) + (t263 * t200 * MDP(19) + (-MDP(20) * t220 + t201 * t267) * t199) * qJD(2);];
tauc = t1;
