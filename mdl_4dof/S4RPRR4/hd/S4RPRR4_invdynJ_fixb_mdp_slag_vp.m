% Calculate vector of inverse dynamics joint torques for
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:39
% EndTime: 2019-12-31 16:50:40
% DurationCPUTime: 1.17s
% Computational Cost: add. (596->195), mult. (1242->287), div. (0->0), fcn. (781->10), ass. (0->91)
t161 = cos(qJ(3));
t189 = qJD(1) * qJD(3);
t184 = t161 * t189;
t158 = sin(qJ(3));
t187 = qJDD(1) * t158;
t213 = qJD(3) * qJD(4) + t184 + t187;
t155 = sin(pkin(7));
t143 = pkin(1) * t155 + pkin(5);
t136 = t143 * qJDD(1);
t138 = t143 * qJD(1);
t125 = qJD(2) * t158 + t138 * t161;
t197 = qJD(3) * t125;
t112 = -qJDD(3) * pkin(3) - qJDD(2) * t161 + t136 * t158 + t197;
t199 = qJD(1) * t161;
t142 = -qJD(4) + t199;
t152 = qJ(1) + pkin(7);
t146 = sin(t152);
t147 = cos(t152);
t176 = g(1) * t147 + g(2) * t146;
t177 = pkin(3) * t158 - pkin(6) * t161;
t212 = (pkin(6) * qJD(4) + t177 * qJD(1)) * t142 + t176 * t158 - g(3) * t161 - t112;
t211 = g(3) * t158;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t178 = t157 * qJDD(3) + t213 * t160;
t192 = qJD(4) * t158;
t183 = qJD(1) * t192;
t114 = -t157 * t183 + t178;
t210 = t114 * t157;
t190 = t160 * qJD(3);
t200 = qJD(1) * t158;
t129 = t157 * t200 - t190;
t209 = t129 * t142;
t195 = qJD(3) * t157;
t131 = t160 * t200 + t195;
t208 = t131 * t142;
t186 = t161 * qJDD(1);
t128 = t158 * t189 + qJDD(4) - t186;
t207 = t157 * t128;
t206 = t157 * t161;
t205 = t160 * t128;
t204 = t160 * t142;
t203 = t160 * t161;
t202 = qJDD(2) - g(3);
t153 = t158 ^ 2;
t201 = -t161 ^ 2 + t153;
t156 = cos(pkin(7));
t144 = -pkin(1) * t156 - pkin(2);
t139 = qJD(1) * t144;
t124 = qJD(2) * t161 - t138 * t158;
t198 = qJD(3) * t124;
t196 = qJD(3) * t129;
t194 = qJD(3) * t158;
t193 = qJD(4) * t157;
t191 = qJD(4) * t160;
t185 = t142 * t195;
t118 = qJD(3) * pkin(6) + t125;
t180 = t142 * t143 + t118;
t111 = qJDD(3) * pkin(6) + qJDD(2) * t158 + t136 * t161 + t198;
t127 = -pkin(3) * t161 - pkin(6) * t158 + t144;
t119 = t127 * qJD(1);
t179 = -qJD(4) * t119 - t111;
t159 = sin(qJ(1));
t162 = cos(qJ(1));
t175 = g(1) * t159 - g(2) * t162;
t174 = qJD(3) * t138 - t202;
t133 = t177 * qJD(3);
t173 = t142 * t191 - t207;
t172 = -t142 * t193 - t205;
t171 = -t157 * t192 + t161 * t190;
t170 = 0.2e1 * t139 * qJD(3) - qJDD(3) * t143;
t117 = -qJD(3) * pkin(3) - t124;
t169 = qJD(3) * t117 - t128 * t143 - t179;
t168 = -pkin(6) * t128 + (-t117 - t124) * t142;
t163 = qJD(3) ^ 2;
t167 = g(1) * t146 - g(2) * t147 - 0.2e1 * qJDD(1) * t144 - t143 * t163;
t166 = -qJD(1) * t139 - qJD(2) * qJD(3) - t136 + t176;
t149 = t160 * qJDD(3);
t135 = qJDD(3) * t161 - t158 * t163;
t134 = qJDD(3) * t158 + t161 * t163;
t126 = t131 * t194;
t123 = t146 * t157 + t147 * t203;
t122 = t146 * t160 - t147 * t206;
t121 = -t146 * t203 + t147 * t157;
t120 = t146 * t206 + t147 * t160;
t116 = qJD(1) * t133 + t127 * qJDD(1);
t115 = t160 * t183 - t149 + (t187 + (qJD(4) + t199) * qJD(3)) * t157;
t113 = t160 * t116;
t110 = t118 * t160 + t119 * t157;
t109 = -t118 * t157 + t119 * t160;
t1 = [qJDD(1) * MDP(1) + t175 * MDP(2) + (g(1) * t162 + g(2) * t159) * MDP(3) + (t175 + (t155 ^ 2 + t156 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t153 + 0.2e1 * t158 * t184) * MDP(5) + 0.2e1 * (t158 * t186 - t201 * t189) * MDP(6) + t134 * MDP(7) + t135 * MDP(8) + (t170 * t158 + t167 * t161) * MDP(10) + (-t167 * t158 + t170 * t161) * MDP(11) + (t114 * t158 * t160 + t171 * t131) * MDP(12) + ((-t129 * t160 - t131 * t157) * t161 * qJD(3) + (-t210 - t115 * t160 + (t129 * t157 - t131 * t160) * qJD(4)) * t158) * MDP(13) + (-t114 * t161 - t171 * t142 + t158 * t205 + t126) * MDP(14) + ((t115 + t185) * t161 + (t173 - t196) * t158) * MDP(15) + (-t128 * t161 - t142 * t194) * MDP(16) + (t158 * t143 * t115 - g(1) * t121 - g(2) * t123 - t113 * t161 + (t161 * t143 * t129 + t109 * t158) * qJD(3) + (t127 * t128 - t133 * t142 + (t117 * t158 + t180 * t161) * qJD(4)) * t160 + (-(-qJD(4) * t127 + t143 * t194) * t142 + t112 * t158 + t169 * t161) * t157) * MDP(17) + ((t127 * t191 + t133 * t157) * t142 - t127 * t207 - g(1) * t120 - g(2) * t122 + (qJD(3) * t143 * t131 + (-t180 * qJD(4) + t116) * t157 + t169 * t160) * t161 + (-t117 * t193 + t112 * t160 + t143 * t114 + (-t143 * t204 - t110) * qJD(3)) * t158) * MDP(18); t202 * MDP(4) + t135 * MDP(10) - t134 * MDP(11) + t126 * MDP(18) + ((-t115 + t185) * MDP(17) + (t142 * t190 - t114) * MDP(18)) * t161 + ((t173 + t196) * MDP(17) + t172 * MDP(18)) * t158; MDP(7) * t187 + MDP(8) * t186 + qJDD(3) * MDP(9) + (t166 * t158 - t174 * t161 + t197) * MDP(10) + (t174 * t158 + t166 * t161 + t198) * MDP(11) + (-t131 * t204 + t210) * MDP(12) + ((t114 + t209) * t160 + (-t115 + t208) * t157) * MDP(13) + ((-t131 * t158 + t142 * t203) * qJD(1) - t173) * MDP(14) + ((t129 * t158 - t142 * t206) * qJD(1) - t172) * MDP(15) + t142 * MDP(16) * t200 + (-pkin(3) * t115 - t109 * t200 - t125 * t129 + t168 * t157 + t212 * t160) * MDP(17) + (-pkin(3) * t114 + t110 * t200 - t125 * t131 - t212 * t157 + t168 * t160) * MDP(18) + (-t158 * t161 * MDP(5) + t201 * MDP(6)) * qJD(1) ^ 2; t131 * t129 * MDP(12) + (-t129 ^ 2 + t131 ^ 2) * MDP(13) + (t178 - t209) * MDP(14) + (t149 - t208) * MDP(15) + t128 * MDP(16) + (-g(1) * t122 + g(2) * t120 - t110 * t142 - t117 * t131 + t113) * MDP(17) + (g(1) * t123 - g(2) * t121 - t109 * t142 + t117 * t129) * MDP(18) + ((-t111 + t211) * MDP(18) + (-MDP(15) * t200 - MDP(17) * t118 - t119 * MDP(18)) * qJD(4)) * t160 + (-MDP(14) * t183 - t213 * MDP(15) + (t179 + t211) * MDP(17) + (qJD(4) * t118 - t116) * MDP(18)) * t157;];
tau = t1;
