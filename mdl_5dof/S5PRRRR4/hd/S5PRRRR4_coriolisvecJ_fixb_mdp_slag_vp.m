% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:08:01
% EndTime: 2019-12-05 17:08:04
% DurationCPUTime: 0.90s
% Computational Cost: add. (689->133), mult. (1229->208), div. (0->0), fcn. (752->6), ass. (0->84)
t243 = pkin(7) + pkin(8);
t182 = sin(qJ(5));
t183 = sin(qJ(4));
t185 = cos(qJ(5));
t186 = cos(qJ(4));
t156 = t182 * t186 + t183 * t185;
t179 = qJD(2) + qJD(3);
t147 = t156 * t179;
t178 = qJD(4) + qJD(5);
t239 = qJD(5) - t178;
t242 = MDP(8) * t183;
t241 = (t183 ^ 2 - t186 ^ 2) * MDP(9);
t184 = sin(qJ(3));
t235 = pkin(2) * qJD(2);
t201 = t243 * t179 + t184 * t235;
t140 = t186 * qJD(1) - t201 * t183;
t240 = MDP(13) * t183 + MDP(14) * t186;
t141 = t183 * qJD(1) + t201 * t186;
t187 = cos(qJ(3));
t237 = pkin(2) * t187;
t172 = pkin(2) * t184 + pkin(7);
t236 = -pkin(8) - t172;
t234 = pkin(2) * qJD(3);
t233 = t178 * t187;
t232 = t179 * t183;
t231 = t182 * t183;
t188 = qJD(4) ^ 2;
t228 = t183 * t188;
t227 = t185 * t141;
t226 = t185 * t186;
t225 = t186 * t188;
t155 = -t226 + t231;
t136 = t178 * t155;
t174 = -pkin(4) * t186 - pkin(3);
t209 = t187 * t235;
t148 = t174 * t179 - t209;
t208 = qJD(2) * t234;
t199 = t184 * t208;
t218 = qJD(4) * t183;
t204 = t179 * t218;
t149 = pkin(4) * t204 + t199;
t224 = -t148 * t136 + t149 * t156;
t137 = t178 * t156;
t223 = t148 * t137 + t149 * t155;
t161 = -pkin(3) * t179 - t209;
t217 = qJD(4) * t186;
t222 = t161 * t217 + t183 * t199;
t216 = qJD(4) * t187;
t214 = -qJD(2) - t179;
t213 = -qJD(3) + t179;
t212 = pkin(4) * t232;
t211 = t187 * t234;
t210 = pkin(4) * t218;
t207 = t179 * t231;
t206 = t179 * t226;
t205 = qJD(4) * t243;
t203 = t179 * t217;
t135 = qJD(4) * pkin(4) + t140;
t202 = -pkin(4) * t178 - t135;
t200 = qJD(4) * t236;
t198 = t187 * t208;
t130 = qJD(4) * t140 + t186 * t198;
t131 = -qJD(4) * t141 - t183 * t198;
t193 = -t182 * t130 + t185 * t131 - t148 * t147;
t128 = qJD(5) * t206 - t178 * t207 + t185 * t203;
t145 = -t206 + t207;
t192 = t147 * t145 * MDP(15) + (-t145 ^ 2 + t147 ^ 2) * MDP(16) + (t145 * t178 + t128) * MDP(17);
t190 = t148 * t145 + (t239 * t141 - t131) * t182;
t129 = t137 * t179;
t189 = -MDP(11) * t228 + (-t128 * t155 - t129 * t156 + t136 * t145 - t137 * t147) * MDP(16) + (t128 * t156 - t136 * t147) * MDP(15) - 0.2e1 * t179 * qJD(4) * t241 + 0.2e1 * t203 * t242 + MDP(10) * t225 + (-t136 * MDP(17) - t137 * MDP(18)) * t178;
t176 = t186 * pkin(8);
t173 = -pkin(3) - t237;
t170 = pkin(7) * t186 + t176;
t169 = t243 * t183;
t164 = t174 - t237;
t159 = t184 * t234 + t210;
t158 = t186 * t205;
t157 = t183 * t205;
t153 = t172 * t186 + t176;
t152 = t236 * t183;
t150 = t161 * t218;
t143 = -t183 * t211 + t186 * t200;
t142 = t183 * t200 + t186 * t211;
t1 = [-t240 * t188 + (-MDP(20) * t137 + MDP(21) * t136) * t178; (-t172 * t225 + t173 * t204 + t150) * MDP(13) + (t172 * t228 + t173 * t203 + t222) * MDP(14) + (t159 * t145 + t164 * t129 + (-t142 * t182 + t143 * t185 + (-t152 * t182 - t153 * t185) * qJD(5)) * t178 + t223) * MDP(20) + (t159 * t147 + t164 * t128 - (t142 * t185 + t143 * t182 + (t152 * t185 - t153 * t182) * qJD(5)) * t178 + t224) * MDP(21) + ((t214 * MDP(7) - qJD(4) * t240) * t187 + (MDP(14) * t232 + (MDP(13) * t186 + MDP(6)) * t214) * t184) * t234 + t189; (-pkin(3) * t204 - pkin(7) * t225 + t150 + (t213 * t186 * t184 + t183 * t216) * t235) * MDP(13) + (-pkin(3) * t203 + pkin(7) * t228 + (-t184 * t232 + t186 * t216) * t235 + t222) * MDP(14) + (t145 * t210 + t174 * t129 + (t157 * t182 - t158 * t185 + (t169 * t182 - t170 * t185) * qJD(5)) * t178 + (-t184 * t145 + t156 * t233) * t235 + t223) * MDP(20) + (t147 * t210 + t174 * t128 - (-t157 * t185 - t158 * t182 + (-t169 * t185 - t170 * t182) * qJD(5)) * t178 + (-t184 * t147 - t155 * t233) * t235 + t224) * MDP(21) + t189 + (t184 * MDP(6) + MDP(7) * t187) * t213 * t235; (-t145 * t212 - (-t140 * t182 - t227) * t178 + (t202 * t182 - t227) * qJD(5) + t193) * MDP(20) + (-t147 * t212 + (t202 * qJD(5) + t140 * t178 - t130) * t185 + t190) * MDP(21) + t192 + t240 * (-t161 * t179 - t198) + (-t186 * t242 + t241) * t179 ^ 2; (t193 + t239 * (-t135 * t182 - t227)) * MDP(20) + ((-t239 * t135 - t130) * t185 + t190) * MDP(21) + t192;];
tauc = t1;
