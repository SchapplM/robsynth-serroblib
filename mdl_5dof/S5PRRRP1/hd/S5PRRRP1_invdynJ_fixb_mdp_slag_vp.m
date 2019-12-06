% Calculate vector of inverse dynamics joint torques for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:20
% EndTime: 2019-12-05 16:40:24
% DurationCPUTime: 0.95s
% Computational Cost: add. (753->170), mult. (1056->224), div. (0->0), fcn. (527->8), ass. (0->91)
t226 = qJ(5) + pkin(7);
t172 = sin(qJ(3));
t217 = pkin(2) * qJD(2);
t200 = t172 * t217;
t174 = cos(qJ(3));
t222 = pkin(2) * t174;
t209 = qJD(3) * t200 - qJDD(2) * t222;
t165 = qJDD(2) + qJDD(3);
t218 = t165 * pkin(3);
t166 = pkin(8) + qJ(2);
t160 = qJ(3) + t166;
t151 = cos(t160);
t220 = g(2) * t151;
t225 = t209 - t218 + t220;
t150 = sin(t160);
t224 = g(1) * t151 + g(2) * t150;
t147 = g(1) * t150;
t223 = t147 - t220;
t167 = qJD(2) + qJD(3);
t221 = pkin(3) * t167;
t173 = cos(qJ(4));
t219 = g(3) * t173;
t216 = qJD(4) * pkin(4);
t171 = sin(qJ(4));
t215 = t165 * t171;
t214 = t171 * t173;
t152 = pkin(2) * t172 + pkin(7);
t213 = -qJ(5) - t152;
t212 = qJDD(1) - g(3);
t192 = t226 * t167 + t200;
t121 = t173 * qJD(1) - t192 * t171;
t117 = t121 + t216;
t211 = t117 - t121;
t168 = t171 ^ 2;
t169 = t173 ^ 2;
t208 = t168 - t169;
t207 = qJD(3) * t172;
t206 = qJD(3) * t174;
t205 = qJD(4) * t171;
t153 = pkin(4) * t173 + pkin(3);
t199 = t174 * t217;
t125 = -t153 * t167 + qJD(5) - t199;
t204 = t125 * MDP(16);
t203 = qJDD(2) * t172;
t202 = qJDD(4) * MDP(13);
t201 = qJDD(4) * MDP(14);
t198 = pkin(2) * t206;
t197 = pkin(4) * t204;
t196 = t167 * t207;
t195 = t167 * t205;
t193 = qJD(4) * t226;
t191 = qJD(4) * t213;
t190 = (-t168 - t169) * MDP(15);
t189 = MDP(16) * t199;
t157 = sin(t166);
t158 = cos(t166);
t187 = g(1) * t157 - g(2) * t158;
t122 = qJD(1) * t171 + t192 * t173;
t186 = t117 * t171 - t122 * t173;
t185 = t192 * qJD(4);
t184 = t199 - t221;
t154 = -pkin(3) - t222;
t183 = t154 * t167 - t198;
t182 = qJD(2) * t206 + t203;
t127 = t182 * pkin(2) + pkin(7) * t165;
t138 = -t199 - t221;
t181 = -t138 * t167 - t127 + t224;
t175 = qJD(4) ^ 2;
t180 = pkin(7) * t175 - t167 * t200 - t218;
t179 = qJ(5) * t165 + qJD(1) * qJD(4) + qJD(5) * t167 + t127;
t116 = pkin(4) * t195 - t153 * t165 + qJDD(5) + t209;
t178 = pkin(2) * t196 + t152 * t175 + t154 * t165;
t177 = -g(1) * (-t150 * t153 + t151 * t226) - g(2) * (t150 * t226 + t151 * t153);
t115 = (qJDD(1) - t185) * t171 + t179 * t173;
t139 = qJDD(4) * t171 + t173 * t175;
t140 = qJDD(4) * t173 - t171 * t175;
t176 = 0.2e1 * (-t208 * t167 * qJD(4) + t165 * t214) * MDP(9) + (t165 * t168 + 0.2e1 * t173 * t195) * MDP(8) + t139 * MDP(10) + t140 * MDP(11) + t165 * MDP(5) + (t138 * t205 + t173 * t147) * MDP(13) + t224 * MDP(7) + (t115 * t173 - t224) * MDP(15) + (t138 * qJD(4) * t173 + t225 * t171) * MDP(14) + (-t209 + t223) * MDP(6);
t164 = t167 ^ 2;
t163 = t173 * qJ(5);
t161 = t173 * qJD(5);
t159 = t173 * qJDD(1);
t144 = pkin(7) * t173 + t163;
t143 = t226 * t171;
t133 = t152 * t173 + t163;
t132 = t213 * t171;
t129 = -qJD(5) * t171 - t173 * t193;
t128 = -t171 * t193 + t161;
t120 = (-qJD(5) - t198) * t171 + t173 * t191;
t119 = t171 * t191 + t173 * t198 + t161;
t114 = qJDD(4) * pkin(4) - t179 * t171 - t173 * t185 + t159;
t1 = [t212 * MDP(1) + t140 * MDP(13) - t139 * MDP(14) + (-t186 * qJD(4) + t114 * t173 + t115 * t171 - g(3)) * MDP(16); t176 + (-t116 * pkin(3) + t114 * t132 + t115 * t133 + t117 * t120 + t122 * t119 + t177) * MDP(16) + ((t165 * t174 - t196) * MDP(6) + (-t165 * t172 - t167 * t206 - t182) * MDP(7) + (-t116 * t174 + t125 * t207 + t187) * MDP(16)) * pkin(2) + t187 * MDP(3) + (g(1) * t158 + g(2) * t157) * MDP(4) + qJDD(2) * MDP(2) + ((-t178 - t225) * MDP(13) - t152 * t201 + (t119 * t167 + t133 * t165) * MDP(15) - t116 * pkin(4) * MDP(16) + (t183 * MDP(14) + (-t132 * t167 - t117) * MDP(15)) * qJD(4)) * t173 + (-t152 * t202 + (t178 - t147) * MDP(14) + (-t120 * t167 - t132 * t165 - t114) * MDP(15) + (t183 * MDP(13) + (-t133 * t167 - t122) * MDP(15) + t197) * qJD(4)) * t171; (-pkin(7) * t202 + (t180 - t147) * MDP(14) + (-t129 * t167 + t143 * t165 - t114) * MDP(15) + t117 * t189 + (t184 * MDP(13) + (-t144 * t167 - t122) * MDP(15) + t197) * qJD(4)) * t171 + ((-t180 - t225) * MDP(13) - pkin(7) * t201 + (t128 * t167 + t144 * t165) * MDP(15) - t122 * t189 + (t184 * MDP(14) + (t143 * t167 - t117) * MDP(15)) * qJD(4)) * t173 + (-MDP(7) * t203 + ((t167 * MDP(6) - t204) * t172 + (-qJD(3) * MDP(7) + (MDP(7) + t190) * t167) * t174) * qJD(2)) * pkin(2) + t176 + (-t114 * t143 + t115 * t144 - t116 * t153 + t117 * t129 + t122 * t128 + t177) * MDP(16); MDP(10) * t215 + t173 * t165 * MDP(11) + qJDD(4) * MDP(12) + (t181 * t171 + t159 - t219) * MDP(13) + (-t212 * t171 + t181 * t173) * MDP(14) + (-pkin(4) * t215 + (t211 - t216) * t173 * t167) * MDP(15) + (t211 * t122 + (-t219 + t114 + (-t125 * t167 + t224) * t171) * pkin(4)) * MDP(16) + (-MDP(8) * t214 + t208 * MDP(9)) * t164; (t186 * t167 + t116 - t223) * MDP(16) + t164 * t190;];
tau = t1;
