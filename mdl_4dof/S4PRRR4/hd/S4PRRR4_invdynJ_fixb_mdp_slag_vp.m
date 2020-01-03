% Calculate vector of inverse dynamics joint torques for
% S4PRRR4
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:42
% EndTime: 2019-12-31 16:32:43
% DurationCPUTime: 0.66s
% Computational Cost: add. (432->132), mult. (889->191), div. (0->0), fcn. (585->8), ass. (0->72)
t148 = pkin(7) + qJ(2);
t141 = sin(t148);
t142 = cos(t148);
t171 = g(1) * t142 + g(2) * t141;
t154 = sin(qJ(3));
t156 = cos(qJ(3));
t193 = pkin(5) + pkin(6);
t177 = qJD(2) * t193;
t110 = t156 * qJD(1) - t154 * t177;
t149 = qJD(3) + qJD(4);
t190 = qJD(3) * pkin(3);
t189 = (qJDD(2) * pkin(2));
t128 = t193 * t156;
t185 = qJD(1) * t154;
t111 = qJD(2) * t128 + t185;
t155 = cos(qJ(4));
t188 = t111 * t155;
t187 = qJDD(1) - g(3);
t150 = t154 ^ 2;
t186 = -t156 ^ 2 + t150;
t184 = qJD(2) * t154;
t183 = qJD(2) * t156;
t153 = sin(qJ(4));
t182 = qJD(4) * t153;
t181 = qJD(2) * qJD(3);
t180 = qJDD(2) * t154;
t179 = qJDD(2) * t156;
t178 = t154 * t190;
t127 = t193 * t154;
t140 = -pkin(3) * t156 - pkin(2);
t176 = qJD(3) * t193;
t175 = t153 * t184;
t174 = t155 * t183;
t173 = t156 * t181;
t170 = g(1) * t141 - g(2) * t142;
t169 = t153 * t180 - t155 * t179;
t116 = t153 * t154 - t155 * t156;
t106 = t149 * t116;
t117 = t153 * t156 + t154 * t155;
t147 = qJDD(3) + qJDD(4);
t168 = -t106 * t149 + t117 * t147;
t109 = t110 + t190;
t167 = -t109 * t153 - t188;
t166 = -t127 * t155 - t128 * t153;
t165 = -t127 * t153 + t128 * t155;
t164 = -0.2e1 * pkin(2) * t181 - pkin(5) * qJDD(3);
t100 = qJD(4) * t174 - t149 * t175 + t153 * t179 + (t173 + t180) * t155;
t113 = -t174 + t175;
t115 = -t153 * t183 - t155 * t184;
t163 = -t115 * t113 * MDP(12) + (-t113 ^ 2 + t115 ^ 2) * MDP(13) + t147 * MDP(16) + (t113 * t149 + t100) * MDP(14) + (-t169 + (-qJD(2) * t117 - t115) * t149) * MDP(15);
t157 = qJD(3) ^ 2;
t162 = -pkin(5) * t157 + t170 + (2 * t189);
t158 = qJD(2) ^ 2;
t161 = pkin(2) * t158 - pkin(5) * qJDD(2) + t171;
t143 = t156 * qJDD(1);
t104 = qJDD(3) * pkin(3) + t143 - qJDD(2) * t127 + (-t156 * t177 - t185) * qJD(3);
t126 = t140 * qJD(2);
t152 = qJ(3) + qJ(4);
t145 = sin(t152);
t146 = cos(t152);
t160 = t126 * t113 + t111 * t182 + g(3) * t145 + (-t111 * t149 - t104) * t153 + t171 * t146;
t107 = t149 * t117;
t105 = t110 * qJD(3) + t154 * qJDD(1) + qJDD(2) * t128;
t159 = -g(3) * t146 + t167 * qJD(4) + t155 * t104 - t153 * t105 + t126 * t115 + t171 * t145;
t125 = qJDD(3) * t156 - t154 * t157;
t124 = qJDD(3) * t154 + t156 * t157;
t119 = t156 * t176;
t118 = t154 * t176;
t112 = -t189 + (t154 * t181 - t179) * pkin(3);
t101 = t107 * qJD(2) + t169;
t99 = -t107 * t149 - t116 * t147;
t1 = [t187 * MDP(1) + t125 * MDP(10) - t124 * MDP(11) + t99 * MDP(17) - t168 * MDP(18); qJDD(2) * MDP(2) + t170 * MDP(3) + t171 * MDP(4) + (qJDD(2) * t150 + 0.2e1 * t154 * t173) * MDP(5) + 0.2e1 * (t154 * t179 - t186 * t181) * MDP(6) + t124 * MDP(7) + t125 * MDP(8) + (t164 * t154 + t162 * t156) * MDP(10) + (-t162 * t154 + t164 * t156) * MDP(11) + (t100 * t117 + t106 * t115) * MDP(12) + (-t100 * t116 - t101 * t117 + t106 * t113 + t107 * t115) * MDP(13) + t168 * MDP(14) + t99 * MDP(15) + (t113 * t178 + t140 * t101 + t112 * t116 + t126 * t107 + (-t165 * qJD(4) + t118 * t153 - t119 * t155) * t149 + t166 * t147 + t170 * t146) * MDP(17) + (-t115 * t178 + t140 * t100 + t112 * t117 - t126 * t106 - (t166 * qJD(4) - t118 * t155 - t119 * t153) * t149 - t165 * t147 - t170 * t145) * MDP(18); MDP(7) * t180 + MDP(8) * t179 + qJDD(3) * MDP(9) + (-g(3) * t156 + t161 * t154 + t143) * MDP(10) + (-t187 * t154 + t161 * t156) * MDP(11) + (-(-t110 * t153 - t188) * t149 + (-t113 * t184 + t147 * t155 - t149 * t182) * pkin(3) + t159) * MDP(17) + ((-qJD(4) * t109 + t110 * t149 - t105) * t155 + (-qJD(4) * t149 * t155 + t115 * t184 - t147 * t153) * pkin(3) + t160) * MDP(18) + t163 + (-t154 * t156 * MDP(5) + t186 * MDP(6)) * t158; (-t167 * t149 + t159) * MDP(17) + ((-t105 + (-qJD(4) + t149) * t109) * t155 + t160) * MDP(18) + t163;];
tau = t1;
