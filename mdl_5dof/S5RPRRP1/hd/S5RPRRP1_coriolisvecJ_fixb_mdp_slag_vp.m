% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:00:01
% EndTime: 2019-12-05 18:00:05
% DurationCPUTime: 0.84s
% Computational Cost: add. (1025->160), mult. (2147->211), div. (0->0), fcn. (1287->4), ass. (0->87)
t198 = cos(qJ(3));
t230 = qJD(3) * t198;
t179 = pkin(3) * t230 + qJD(2);
t190 = qJD(3) + qJD(4);
t196 = sin(qJ(3));
t252 = (t198 * MDP(12) - t196 * MDP(13)) * qJ(2) - t198 * t196 * MDP(7) + (t196 ^ 2 - t198 ^ 2) * MDP(8);
t199 = -pkin(1) - pkin(6);
t178 = t199 * qJD(1) + qJD(2);
t232 = qJD(1) * t198;
t156 = -pkin(7) * t232 + t198 * t178;
t152 = qJD(3) * pkin(3) + t156;
t219 = pkin(7) * qJD(1) - t178;
t154 = t219 * t230;
t197 = cos(qJ(4));
t251 = (qJD(4) * t152 - t154) * t197;
t195 = sin(qJ(4));
t233 = qJD(1) * t196;
t221 = t195 * t233;
t222 = t197 * t232;
t164 = -t221 + t222;
t157 = t164 * qJ(5);
t155 = -pkin(7) * t233 + t178 * t196;
t149 = t195 * t155;
t217 = t197 * t152 - t149;
t250 = t157 - t217;
t249 = t196 * MDP(12) + t198 * MDP(13);
t169 = t195 * t198 + t196 * t197;
t162 = t169 * qJD(1);
t246 = t164 ^ 2;
t245 = pkin(7) - t199;
t243 = qJ(5) * t162;
t175 = pkin(3) * t233 + qJD(1) * qJ(2);
t145 = pkin(4) * t162 + qJD(5) + t175;
t242 = t145 * t164;
t174 = t245 * t198;
t241 = t174 * t197;
t240 = t175 * t164;
t150 = t197 * t155;
t239 = t197 * t198;
t184 = t196 * pkin(3) + qJ(2);
t131 = pkin(4) * t190 - t250;
t238 = t131 + t250;
t237 = t197 * t156 - t149;
t236 = t190 * t221;
t172 = t179 * qJD(1);
t231 = qJD(3) * t196;
t229 = qJD(4) * t195;
t224 = pkin(3) * t232;
t220 = -pkin(3) * t190 - t152;
t218 = -qJ(2) * MDP(6) - MDP(5);
t153 = t219 * t231;
t216 = t197 * t153 + t195 * t154;
t215 = t195 * t153 - t155 * t229;
t214 = -t156 * t195 - t150;
t211 = t190 * t239;
t142 = qJD(1) * t211 - t236;
t212 = pkin(4) * t142 + t172;
t143 = t190 * t169;
t141 = t143 * qJD(1);
t170 = -t195 * t196 + t239;
t210 = -t141 * t170 - t143 * t164;
t209 = -t152 * t195 - t150;
t173 = t245 * t196;
t208 = t173 * t197 + t174 * t195;
t207 = t175 * t162 - t215;
t161 = t162 ^ 2;
t206 = t164 * t162 * MDP(14) + (t236 + (t164 - t222) * t190) * MDP(17) + (-t161 + t246) * MDP(15);
t167 = t245 * t231;
t168 = qJD(3) * t174;
t205 = qJD(4) * t241 - t195 * t167 + t197 * t168 - t173 * t229;
t127 = -qJ(5) * t142 - qJD(5) * t162 + t215 + t251;
t203 = t209 * qJD(4) + t216;
t128 = qJ(5) * t141 - qJD(5) * t164 + t203;
t133 = -t209 - t243;
t144 = -t195 * t231 - t196 * t229 + t211;
t204 = t127 * t169 + t128 * t170 - t131 * t143 + t133 * t144;
t202 = t208 * qJD(4) + t197 * t167 + t168 * t195;
t201 = qJD(1) ^ 2;
t200 = qJD(3) ^ 2;
t185 = pkin(3) * t197 + pkin(4);
t140 = -qJ(5) * t169 - t208;
t139 = -qJ(5) * t170 + t173 * t195 - t241;
t135 = -t157 + t237;
t134 = t214 + t243;
t130 = qJ(5) * t143 - qJD(5) * t170 + t202;
t129 = -qJ(5) * t144 - qJD(5) * t169 - t205;
t1 = [t210 * MDP(14) + (t141 * t169 - t142 * t170 + t143 * t162 - t144 * t164) * MDP(15) + (t184 * t142 + t175 * t144 + t179 * t162 + t172 * t169) * MDP(19) + (-t184 * t141 - t175 * t143 + t179 * t164 + t172 * t170) * MDP(20) + (-t129 * t162 - t130 * t164 + t139 * t141 - t140 * t142 - t204) * MDP(21) + (t127 * t140 + t133 * t129 + t128 * t139 + t131 * t130 + t212 * (pkin(4) * t169 + t184) + t145 * (pkin(4) * t144 + t179)) * MDP(22) + (-t143 * MDP(16) - t144 * MDP(17) + t202 * MDP(19) + t205 * MDP(20)) * t190 + ((-MDP(13) * t199 - MDP(10)) * t198 + (-MDP(12) * t199 - MDP(9)) * t196) * t200 + (0.2e1 * (-t218 + t249) * qJD(2) + 0.2e1 * t252 * qJD(3)) * qJD(1); (-qJD(1) * t162 - t143 * t190) * MDP(19) + (-qJD(1) * t164 - t144 * t190) * MDP(20) + (-t142 * t169 - t144 * t162 - t210) * MDP(21) + (-qJD(1) * t145 + t204) * MDP(22) + t218 * t201 + t249 * (-t200 - t201); (-t162 * t224 - t240 - t214 * t190 + (t220 * t195 - t150) * qJD(4) + t216) * MDP(19) + (-t164 * t224 + t237 * t190 + (t220 * qJD(4) + t154) * t197 + t207) * MDP(20) + (t141 * t185 + (t133 + t134) * t164 + (-t131 + t135) * t162 + (-t142 * t195 + (-t162 * t197 + t164 * t195) * qJD(4)) * pkin(3)) * MDP(21) + (-pkin(4) * t242 + t128 * t185 - t131 * t134 - t133 * t135 + (-t145 * t232 + t127 * t195 + (-t131 * t195 + t133 * t197) * qJD(4)) * pkin(3)) * MDP(22) + t206 - t252 * t201; (-t209 * t190 + t203 - t240) * MDP(19) + (t217 * t190 + t207 - t251) * MDP(20) + (pkin(4) * t141 - t238 * t162) * MDP(21) + (t238 * t133 + (t128 - t242) * pkin(4)) * MDP(22) + t206; (-t161 - t246) * MDP(21) + (t131 * t164 + t133 * t162 + t212) * MDP(22);];
tauc = t1;
