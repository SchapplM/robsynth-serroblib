% Calculate vector of inverse dynamics joint torques for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:40
% EndTime: 2019-12-31 17:38:42
% DurationCPUTime: 0.81s
% Computational Cost: add. (481->154), mult. (852->203), div. (0->0), fcn. (582->8), ass. (0->72)
t178 = cos(qJ(2));
t165 = t178 * qJDD(1);
t176 = sin(qJ(2));
t209 = qJD(2) * t176;
t202 = qJD(1) * t209 + qJDD(3) - t165;
t221 = pkin(2) + pkin(3);
t127 = -t221 * qJDD(2) + t202;
t163 = t176 * qJDD(1);
t206 = qJDD(2) * qJ(3);
t208 = t178 * qJD(1);
t132 = t206 + t163 + (qJD(3) + t208) * qJD(2);
t171 = sin(pkin(8));
t173 = cos(pkin(8));
t118 = t127 * t173 - t171 * t132;
t210 = qJD(1) * t176;
t135 = t171 * t208 - t173 * t210;
t192 = t176 * t171 + t178 * t173;
t216 = t171 * t178;
t193 = -t173 * t176 + t216;
t172 = sin(pkin(7));
t174 = cos(pkin(7));
t195 = g(1) * t174 + g(2) * t172;
t226 = (qJD(3) * t171 - t135) * qJD(2) - g(3) * t192 - t195 * t193 - t118;
t119 = t171 * t127 + t173 * t132;
t225 = -g(3) * t193 + t195 * t192 - t119;
t145 = -qJ(3) * t171 - t173 * t221;
t142 = pkin(4) - t145;
t146 = t173 * qJ(3) - t171 * t221;
t143 = -pkin(6) + t146;
t180 = qJD(5) ^ 2;
t223 = -t143 * t180 + (pkin(4) + t142) * qJDD(2) + t226;
t222 = g(3) * t176 + t195 * t178 - t163;
t218 = qJ(3) * t178;
t217 = qJDD(2) * pkin(2);
t215 = t172 * t176;
t214 = t174 * t176;
t177 = cos(qJ(5));
t181 = qJD(2) ^ 2;
t213 = t177 * t181;
t212 = t178 * pkin(2) + t176 * qJ(3);
t175 = sin(qJ(5));
t169 = t175 ^ 2;
t211 = -t177 ^ 2 + t169;
t207 = qJD(2) * qJD(5);
t205 = qJDD(2) * t177;
t204 = 0.2e1 * t207;
t203 = -g(1) * t214 - g(2) * t215 + g(3) * t178;
t136 = t192 * qJD(1);
t199 = qJD(3) * t173 - t136;
t198 = t177 * t204;
t197 = qJD(3) - t208;
t196 = t165 - t203;
t141 = -t221 * qJD(2) + t197;
t152 = qJD(2) * qJ(3) + t210;
t122 = t141 * t173 - t152 * t171;
t194 = (-qJD(2) * pkin(2) + t197) * t176 + t152 * t178;
t148 = qJDD(5) * t177 - t175 * t180;
t147 = -qJDD(5) * t175 - t177 * t180;
t137 = -qJD(2) * t216 + t173 * t209;
t190 = -qJD(2) * t137 + qJDD(2) * t192;
t189 = g(1) * t172 - g(2) * t174 + qJDD(4);
t133 = t202 - t217;
t186 = -t180 * t193 - t190;
t138 = t192 * qJD(2);
t185 = -qJD(5) * t138 + qJDD(5) * t193 - t192 * t207;
t120 = qJD(2) * pkin(4) - t122;
t184 = qJDD(2) * pkin(6) + qJD(2) * t120 + t225;
t183 = -qJDD(5) * t143 + (-qJD(2) * t142 - t120 - t199) * qJD(5);
t156 = t174 * t218;
t155 = t172 * t218;
t123 = t171 * t141 + t173 * t152;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t194 * qJD(2) + t132 * t176 - t133 * t178 - g(3)) * MDP(7) + t190 * MDP(8) + (qJD(2) * t138 - qJDD(2) * t193) * MDP(9) + (-t118 * t192 - t119 * t193 + t122 * t137 + t123 * t138 - g(3)) * MDP(10) + (MDP(3) + MDP(5)) * (qJDD(2) * t178 - t176 * t181) + (-MDP(4) + MDP(6)) * (qJDD(2) * t176 + t178 * t181) + (-t186 * MDP(16) + t185 * MDP(17)) * t177 + (t185 * MDP(16) + t186 * MDP(17)) * t175; qJDD(2) * MDP(2) + t196 * MDP(3) + t222 * MDP(4) + (-qJDD(3) + t196 + 0.2e1 * t217) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t206 - t222) * MDP(6) + (t132 * qJ(3) + t152 * qJD(3) - t133 * pkin(2) - g(1) * (-pkin(2) * t214 + t156) - g(2) * (-pkin(2) * t215 + t155) - g(3) * t212 - t194 * qJD(1)) * MDP(7) + (-qJDD(2) * t145 + t226) * MDP(8) + (t199 * qJD(2) + qJDD(2) * t146 - t225) * MDP(9) + (t119 * t146 + t118 * t145 - t123 * t136 + t122 * t135 - g(1) * t156 - g(2) * t155 - g(3) * (pkin(3) * t178 + t212) + (-t122 * t171 + t123 * t173) * qJD(3) + t195 * t176 * t221) * MDP(10) + (qJDD(2) * t169 + t175 * t198) * MDP(11) + 0.2e1 * (t175 * t205 - t211 * t207) * MDP(12) + t147 * MDP(13) - t148 * MDP(14) + (t183 * t175 + t223 * t177) * MDP(16) + (-t223 * t175 + t183 * t177) * MDP(17); -qJDD(2) * MDP(5) - t181 * MDP(6) + (-qJD(2) * t152 + t133 + t203) * MDP(7) + t203 * MDP(10) + (-qJDD(2) * MDP(8) - t181 * MDP(9) + (-qJD(2) * t123 + t118) * MDP(10) + (t175 * t204 - t205) * MDP(16) + (qJDD(2) * t175 + t198) * MDP(17)) * t173 + (-t181 * MDP(8) + qJDD(2) * MDP(9) + (qJD(2) * t122 + t119) * MDP(10) + (t147 - t213) * MDP(16) + (t175 * t181 - t148) * MDP(17)) * t171; t189 * MDP(10) + t148 * MDP(16) + t147 * MDP(17); qJDD(5) * MDP(15) + t211 * MDP(12) * t181 + (-qJDD(2) * MDP(14) + t189 * MDP(16) + t184 * MDP(17)) * t177 + (-MDP(11) * t213 - qJDD(2) * MDP(13) + t184 * MDP(16) - t189 * MDP(17)) * t175;];
tau = t1;
