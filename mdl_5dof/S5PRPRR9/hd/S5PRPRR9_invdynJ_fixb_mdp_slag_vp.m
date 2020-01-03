% Calculate vector of inverse dynamics joint torques for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:48
% EndTime: 2019-12-31 17:39:49
% DurationCPUTime: 0.50s
% Computational Cost: add. (488->108), mult. (648->137), div. (0->0), fcn. (344->6), ass. (0->56)
t163 = pkin(8) + qJ(2);
t126 = cos(t163);
t160 = sin(t163);
t191 = g(1) * t160 - g(2) * t126;
t164 = qJD(2) - qJD(4);
t190 = t164 ^ 2;
t139 = -pkin(2) - pkin(3);
t171 = qJ(3) * qJD(2);
t189 = -qJD(4) * t171 + t139 * qJDD(2) + qJDD(3);
t188 = -qJDD(3) + t191;
t169 = qJD(5) * t164;
t118 = t139 * qJD(2) + qJD(3);
t165 = (qJD(2) * qJD(3));
t166 = (qJ(3) * qJDD(2));
t187 = qJD(4) * t118 + t165 + t166;
t179 = pkin(2) * qJDD(2);
t185 = t179 + t188;
t136 = sin(qJ(4));
t138 = cos(qJ(4));
t173 = t138 * qJ(3) + t136 * t139;
t105 = -t126 * t138 - t160 * t136;
t106 = t126 * t136 - t160 * t138;
t142 = g(1) * t105 + g(2) * t106 + t189 * t136 + t187 * t138;
t146 = -g(1) * t106 + g(2) * t105 + t187 * t136 - t189 * t138;
t103 = t138 * t118 - t136 * t171;
t180 = t164 * pkin(4);
t98 = -t103 + t180;
t184 = (t103 + t98 + t180) * qJD(5) - pkin(7) * qJDD(5);
t135 = sin(qJ(5));
t137 = cos(qJ(5));
t183 = t137 * MDP(16) - t135 * MDP(17) + MDP(9);
t130 = qJDD(2) - qJDD(4);
t181 = t130 * pkin(4);
t177 = (t136 * qJD(3) + t173 * qJD(4)) * t164;
t176 = (t136 * t118 + t138 * t171) * t164;
t175 = t135 * t137;
t174 = qJDD(1) - g(3);
t132 = t135 ^ 2;
t172 = -t137 ^ 2 + t132;
t156 = -t136 * qJ(3) + t138 * t139;
t140 = qJD(5) ^ 2;
t113 = qJDD(5) * t137 - t140 * t135;
t112 = qJDD(5) * t135 + t140 * t137;
t155 = -t181 - t146;
t151 = g(1) * t126 + g(2) * t160;
t150 = t130 * pkin(7) + t164 * t98 - t142;
t101 = t138 * qJD(3) + t156 * qJD(4);
t109 = pkin(4) - t156;
t110 = -pkin(7) + t173;
t148 = -qJDD(5) * t110 + (-t109 * t164 - t101 - t98) * qJD(5);
t147 = -t151 + (2 * t165);
t145 = pkin(7) * t140 - t155 + t176 + t181;
t144 = -t109 * t130 + t110 * t140 + t155 - t177;
t143 = (-t132 * t130 - 0.2e1 * t169 * t175) * MDP(11) + 0.2e1 * (-t130 * t175 + t172 * t169) * MDP(12) + t112 * MDP(13) + t113 * MDP(14) - t130 * MDP(8);
t141 = qJD(2) ^ 2;
t1 = [-t113 * MDP(16) + t112 * MDP(17) + (MDP(1) + MDP(7)) * t174; qJDD(2) * MDP(2) + t191 * MDP(3) + t151 * MDP(4) + (0.2e1 * t179 + t188) * MDP(5) + (t147 + (2 * t166)) * MDP(6) + (t185 * pkin(2) + (t147 + t166) * qJ(3)) * MDP(7) + (-t156 * t130 + t146 + t177) * MDP(9) + (t101 * t164 + t173 * t130 + t142) * MDP(10) + (t148 * t135 - t144 * t137) * MDP(16) + (t144 * t135 + t148 * t137) * MDP(17) - t143; -qJDD(2) * MDP(5) - t141 * MDP(6) + (-t141 * qJ(3) - t185) * MDP(7) + (-MDP(10) * t190 - t183 * t130 + 0.2e1 * (MDP(16) * t135 + MDP(17) * t137) * t169) * t138 + (t130 * MDP(10) - t112 * MDP(16) - t113 * MDP(17) - t183 * t190) * t136; (-t146 - t176) * MDP(9) + (-t103 * t164 - t142) * MDP(10) + (-t145 * MDP(16) + t184 * MDP(17)) * t137 + (t184 * MDP(16) + t145 * MDP(17)) * t135 + t143; qJDD(5) * MDP(15) + t172 * MDP(12) * t190 + (-t130 * MDP(14) - t174 * MDP(16) + t150 * MDP(17)) * t137 + (-MDP(11) * t137 * t190 - t130 * MDP(13) + t150 * MDP(16) + t174 * MDP(17)) * t135;];
tau = t1;
