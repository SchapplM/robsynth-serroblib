% Calculate vector of inverse dynamics joint torques for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:59
% EndTime: 2019-12-31 17:37:01
% DurationCPUTime: 0.82s
% Computational Cost: add. (436->172), mult. (821->227), div. (0->0), fcn. (569->6), ass. (0->86)
t176 = sin(pkin(8));
t173 = t176 ^ 2;
t177 = cos(pkin(8));
t224 = t177 ^ 2 + t173;
t175 = pkin(7) + qJ(2);
t172 = cos(t175);
t165 = g(2) * t172;
t171 = sin(t175);
t220 = g(1) * t171;
t191 = -t165 + t220;
t222 = MDP(8) + MDP(12);
t221 = pkin(3) * t177;
t219 = -pkin(6) + qJ(3);
t200 = qJDD(2) * t177;
t204 = qJD(2) * qJD(3);
t130 = qJ(3) * t200 + t176 * qJDD(1) + t177 * t204;
t218 = t130 * t176 - g(3);
t217 = qJDD(2) * pkin(2);
t126 = t130 * t177;
t210 = qJ(3) * qJD(2);
t216 = (t176 * qJD(1) + t177 * t210) * t177;
t215 = t171 * t177;
t179 = cos(qJ(5));
t214 = t176 * t179;
t178 = sin(qJ(5));
t213 = t177 * t178;
t212 = qJ(3) * t126 + qJD(3) * t216;
t211 = t172 * pkin(2) + t171 * qJ(3);
t209 = qJD(2) * t178;
t208 = qJD(2) * t179;
t207 = qJD(4) * t176;
t206 = t179 * MDP(18);
t203 = qJD(2) * qJD(4);
t193 = qJ(4) * t176 + pkin(2);
t145 = -t193 - t221;
t202 = qJDD(2) * t145;
t201 = qJDD(2) * t176;
t199 = qJDD(2) * t178;
t198 = qJDD(2) * t179;
t197 = t176 * t208;
t196 = qJD(5) * t214;
t195 = t177 * t209;
t194 = t176 * t203;
t192 = -pkin(2) * t171 + t172 * qJ(3);
t143 = qJD(1) * t177 - t176 * t210;
t190 = g(1) * t172 + g(2) * t171;
t151 = t176 * t198;
t188 = -t177 * t199 + t151;
t129 = -qJ(3) * t201 + qJDD(1) * t177 - t176 * t204;
t187 = qJDD(3) - t194;
t147 = t219 * t176;
t148 = t219 * t177;
t186 = t147 * t179 - t148 * t178;
t185 = t147 * t178 + t148 * t179;
t142 = -t213 + t214;
t141 = t176 * t178 + t177 * t179;
t168 = qJDD(3) - t217;
t184 = -t168 + t217 - t165;
t118 = t187 + t202;
t183 = -t118 - t202 - t165;
t127 = qJDD(4) - t129;
t134 = t141 * qJD(5);
t139 = (pkin(3) + pkin(4)) * t177 + t193;
t182 = t126 - t190 + t224 * (qJ(3) * qJDD(2) + t204);
t181 = qJD(2) ^ 2;
t180 = qJD(5) ^ 2;
t150 = qJD(5) * t195;
t149 = g(1) * t215;
t140 = qJD(4) - t143;
t138 = -t195 + t197;
t136 = t141 * qJD(2);
t135 = -qJD(5) * t213 + t196;
t132 = t145 * qJD(2) + qJD(3);
t124 = t141 * t172;
t123 = t142 * t172;
t122 = t141 * t171;
t121 = t142 * t171;
t120 = t139 * qJD(2) - qJD(3);
t117 = -pkin(6) * t200 + t130;
t116 = -pkin(6) * t201 + t127;
t115 = t139 * qJDD(2) - t187;
t114 = -qJD(5) * t135 - qJDD(5) * t141;
t113 = -qJD(5) * t134 + qJDD(5) * t142;
t112 = qJD(2) * t196 + t141 * qJDD(2) - t150;
t111 = -qJD(2) * t134 + t188;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t129 * t177 + t218) * MDP(8) + (-t127 * t177 + t218) * MDP(12) + t114 * MDP(18) - t113 * MDP(19); qJDD(2) * MDP(2) + t191 * MDP(3) + t190 * MDP(4) + (t184 * t177 + t149) * MDP(5) + (-t184 - t220) * t176 * MDP(6) + (-t129 * t176 + t182) * MDP(7) + (-t168 * pkin(2) - g(1) * t192 - g(2) * t211 + (-t129 * qJ(3) - t143 * qJD(3)) * t176 + t212) * MDP(8) + (t149 + (t183 + t194) * t177) * MDP(9) + (t127 * t176 + t182) * MDP(10) + (t173 * t203 + (t183 + t220) * t176) * MDP(11) + (t118 * t145 - g(1) * (-pkin(3) * t215 + t192) - g(2) * (t172 * t221 + t211) + (t127 * qJ(3) + t191 * qJ(4) + t140 * qJD(3) - t132 * qJD(4)) * t176 + t212) * MDP(12) + (t111 * t142 - t134 * t138) * MDP(13) + (-t111 * t141 - t112 * t142 + t134 * t136 - t135 * t138) * MDP(14) + t113 * MDP(15) + t114 * MDP(16) + (t136 * t207 + t139 * t112 + t115 * t141 + t120 * t135 + t186 * qJDD(5) + g(1) * t122 - g(2) * t124 + (t142 * qJD(3) - t185 * qJD(5)) * qJD(5)) * MDP(18) + (t138 * t207 + t139 * t111 + t115 * t142 - t120 * t134 - t185 * qJDD(5) + g(1) * t121 - g(2) * t123 + (-t141 * qJD(3) - t186 * qJD(5)) * qJD(5)) * MDP(19); (-qJD(5) * t138 + t150) * MDP(18) + (qJD(5) * t136 - t151) * MDP(19) - (MDP(7) + MDP(10)) * t224 * t181 + ((t143 * t176 - t216) * MDP(8) + (-t140 * t176 - t207 - t216) * MDP(12) + (t141 * MDP(19) - t176 * t206) * qJD(5)) * qJD(2) + (-t222 * pkin(2) + (-qJ(4) * MDP(12) - t178 * MDP(18) - MDP(11) + MDP(6)) * t176 + (-pkin(3) * MDP(12) + t178 * MDP(19) - MDP(5) - MDP(9) - t206) * t177) * qJDD(2) + t222 * (qJDD(3) - t191); -t173 * t181 * MDP(11) + (g(3) * t177 + t127) * MDP(12) + (t179 * qJDD(5) - t178 * t180) * MDP(18) + (-t178 * qJDD(5) - t179 * t180) * MDP(19) + (-t181 * t177 * MDP(9) + qJDD(2) * MDP(10) - t190 * MDP(12) + (t132 * MDP(12) - t136 * MDP(18) - t138 * MDP(19)) * qJD(2)) * t176; t138 * t136 * MDP(13) + (-t136 ^ 2 + t138 ^ 2) * MDP(14) + t188 * MDP(15) + (-t176 * t199 - t177 * t198 + t150) * MDP(16) + qJDD(5) * MDP(17) + (-g(1) * t123 - g(2) * t121 + g(3) * t141 + t179 * t116 - t178 * t117 - t120 * t138) * MDP(18) + (g(1) * t124 + g(2) * t122 + g(3) * t142 - t178 * t116 - t179 * t117 + t120 * t136) * MDP(19) + ((-t176 * t209 - t177 * t208 + t136) * MDP(15) + (t138 - t197) * MDP(16)) * qJD(5);];
tau = t1;
