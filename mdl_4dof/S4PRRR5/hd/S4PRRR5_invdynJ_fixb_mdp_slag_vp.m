% Calculate vector of inverse dynamics joint torques for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:48
% EndTime: 2019-12-31 16:33:49
% DurationCPUTime: 0.58s
% Computational Cost: add. (380->112), mult. (660->162), div. (0->0), fcn. (458->10), ass. (0->61)
t141 = sin(pkin(7));
t142 = cos(pkin(7));
t165 = g(1) * t142 + g(2) * t141;
t148 = cos(qJ(2));
t132 = t148 * qJDD(1);
t145 = sin(qJ(2));
t171 = qJD(1) * qJD(2);
t114 = qJDD(2) * pkin(2) - t145 * t171 + t132;
t140 = qJ(2) + qJ(3);
t134 = cos(t140);
t147 = cos(qJ(3));
t192 = -g(3) * t134 + t147 * t114;
t133 = sin(t140);
t144 = sin(qJ(3));
t174 = qJD(1) * t145;
t168 = qJD(3) * t174;
t166 = g(3) * t133 + t165 * t134 + t144 * t168;
t127 = qJD(2) * pkin(2) + qJD(1) * t148;
t155 = qJD(3) * t127 + qJDD(1) * t145 + t148 * t171;
t136 = qJDD(2) + qJDD(3);
t184 = pkin(3) * t136;
t190 = t155 * t144 + t147 * t168 - t184 - t192;
t189 = t165 * t133;
t137 = qJD(2) + qJD(3);
t188 = t137 ^ 2;
t115 = t144 * t145 - t147 * t148;
t113 = t115 * qJD(1);
t129 = pkin(2) * t144 + pkin(6);
t130 = -pkin(2) * t147 - pkin(3);
t180 = pkin(2) * qJD(3);
t186 = -qJDD(4) * t129 + (t130 * t137 - t147 * t180 - t113) * qJD(4);
t185 = pkin(2) * t136;
t179 = (t127 * t144 + t147 * t174) * t137;
t116 = t144 * t148 + t145 * t147;
t178 = t116 * qJD(1) * t137;
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t177 = t143 * t146;
t176 = qJDD(1) - g(3);
t138 = t143 ^ 2;
t175 = -t146 ^ 2 + t138;
t173 = qJD(4) * t137;
t172 = qJD(4) * t146;
t169 = t137 * t180;
t164 = g(1) * t141 - g(2) * t142;
t163 = -t115 * t136 - t188 * t116;
t162 = t189 + t192;
t109 = t127 * t147 - t144 * t174;
t149 = qJD(4) ^ 2;
t160 = -pkin(6) * t149 + t179 + t184;
t159 = t116 * t149 - t163;
t157 = -pkin(3) * t173 - pkin(6) * qJDD(4) + qJD(4) * t109;
t103 = t115 * t137;
t156 = qJD(4) * t103 - qJDD(4) * t116 + t115 * t173;
t154 = -t129 * t149 - t130 * t136 - t144 * t169 + t178;
t107 = -pkin(3) * t137 - t109;
t153 = -pkin(6) * t136 - t107 * t137 - t114 * t144 - t155 * t147 + t166;
t152 = -t155 - t169;
t151 = 0.2e1 * (t136 * t177 - t175 * t173) * MDP(9) + (0.2e1 * t137 * t143 * t172 + t136 * t138) * MDP(8) + (qJDD(4) * t143 + t146 * t149) * MDP(10) + (qJDD(4) * t146 - t143 * t149) * MDP(11) + t136 * MDP(5) + (t107 * qJD(4) * t143 + t189 * t146) * MDP(13) + (t107 * t172 + t190 * t143) * MDP(14);
t150 = qJD(2) ^ 2;
t1 = [t176 * MDP(1) + (qJDD(2) * t148 - t145 * t150) * MDP(3) + (-qJDD(2) * t145 - t148 * t150) * MDP(4) + t163 * MDP(6) + (t103 * t137 - t116 * t136) * MDP(7) + (-t159 * MDP(13) + t156 * MDP(14)) * t146 + (t156 * MDP(13) + t159 * MDP(14)) * t143; ((-t168 + t185) * MDP(6) + t152 * MDP(7)) * t147 + (t152 * MDP(6) + (-t114 - t185) * MDP(7)) * t144 + ((-t154 - t189) * MDP(14) + t186 * MDP(13)) * t143 + ((t154 - t190) * MDP(13) + t186 * MDP(14)) * t146 + t132 * MDP(3) + (-g(3) * MDP(3) + t165 * MDP(4)) * t148 + (t165 * MDP(3) - t176 * MDP(4)) * t145 + (t162 + t178) * MDP(6) + (-t113 * t137 + t166) * MDP(7) + t151 + qJDD(2) * MDP(2); (t162 + t179) * MDP(6) + (t109 * t137 + t166) * MDP(7) + (-MDP(6) * t168 - t155 * MDP(7)) * t147 + (-t155 * MDP(6) - t114 * MDP(7)) * t144 + ((t160 - t190) * MDP(13) + t157 * MDP(14)) * t146 + (t157 * MDP(13) + (-t160 - t189) * MDP(14)) * t143 + t151; qJDD(4) * MDP(12) + (t153 * t143 - t164 * t146) * MDP(13) + (t164 * t143 + t153 * t146) * MDP(14) + (t143 * MDP(10) + t146 * MDP(11)) * t136 + (-MDP(8) * t177 + t175 * MDP(9)) * t188;];
tau = t1;
