% Calculate vector of inverse dynamics joint torques for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:18
% EndTime: 2019-12-31 17:45:20
% DurationCPUTime: 0.97s
% Computational Cost: add. (482->154), mult. (786->192), div. (0->0), fcn. (498->12), ass. (0->78)
t168 = sin(pkin(7));
t146 = pkin(1) * t168 + qJ(3);
t205 = qJDD(1) * t146;
t167 = sin(pkin(8));
t169 = cos(pkin(8));
t217 = t167 * MDP(8) + t169 * MDP(9);
t204 = t167 ^ 2 + t169 ^ 2;
t162 = qJ(1) + pkin(7);
t155 = sin(t162);
t157 = cos(t162);
t216 = -g(1) * t155 + g(2) * t157;
t170 = cos(pkin(7));
t149 = -pkin(1) * t170 - pkin(2);
t215 = qJDD(1) * t149;
t143 = -qJ(4) + t149;
t214 = qJD(1) * qJD(4) - qJDD(1) * t143;
t212 = pkin(4) * t167;
t210 = -pkin(6) + t143;
t171 = sin(qJ(5));
t208 = t167 * t171;
t173 = cos(qJ(5));
t206 = t169 * t173;
t138 = t146 * qJD(1);
t121 = qJDD(3) - t214;
t115 = t169 * qJDD(2) + t167 * t121;
t203 = qJD(1) * t171;
t202 = qJD(1) * t173;
t164 = qJD(3) * qJD(1);
t200 = qJDD(1) * t167;
t199 = qJDD(1) * t171;
t198 = qJDD(1) * t173;
t174 = cos(qJ(1));
t197 = t174 * pkin(1) + t157 * pkin(2) + t155 * qJ(3);
t196 = -t164 - t205;
t195 = t169 * t202;
t194 = qJD(5) * t206;
t193 = t167 * t203;
t136 = qJD(4) + t138;
t172 = sin(qJ(1));
t192 = -pkin(1) * t172 + t157 * qJ(3);
t190 = t204 * MDP(10);
t114 = -qJDD(2) * t167 + t169 * t121;
t132 = qJDD(4) - t196;
t189 = -g(1) * t157 - g(2) * t155;
t188 = g(1) * t172 - g(2) * t174;
t187 = qJDD(3) + t216;
t140 = t169 * t198;
t186 = -t167 * t199 + t140;
t185 = t114 * t169 + t115 * t167;
t184 = t204 * (qJD(1) * t143 + qJD(3));
t125 = t210 * t167;
t126 = t210 * t169;
t183 = t125 * t173 + t126 * t171;
t182 = t125 * t171 - t126 * t173;
t134 = t167 * t173 + t169 * t171;
t181 = -t206 + t208;
t130 = t134 * qJD(5);
t108 = -qJD(5) * t130 - qJDD(5) * t181;
t131 = -qJD(5) * t208 + t194;
t109 = -qJD(5) * t131 - qJDD(5) * t134;
t180 = -t167 * t202 - t169 * t203;
t179 = t189 + t205;
t178 = -t185 - t216;
t175 = qJD(1) ^ 2;
t161 = pkin(8) + qJ(5);
t156 = cos(t161);
t154 = sin(t161);
t139 = qJD(5) * t193;
t137 = t146 + t212;
t129 = -t193 + t195;
t127 = t134 * qJD(1);
t124 = qJD(1) * t212 + t136;
t120 = pkin(4) * t200 + t132;
t111 = -pkin(6) * t200 + t115;
t110 = -pkin(6) * qJDD(1) * t169 + t114;
t107 = qJD(1) * t194 + qJDD(1) * t134 - t139;
t106 = -qJD(1) * t130 + t186;
t1 = [qJDD(1) * MDP(1) + t188 * MDP(2) + (g(1) * t174 + g(2) * t172) * MDP(3) + (t188 + (t168 ^ 2 + t170 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t187 + 0.2e1 * t215) * MDP(5) + (0.2e1 * t164 + t179 + t205) * MDP(6) + (-t196 * t146 + t138 * qJD(3) + (qJDD(3) + t215) * t149 - g(1) * (-pkin(2) * t155 + t192) - g(2) * t197) * MDP(7) + (t214 * t204 + t178) * MDP(10) + (t132 * t146 + t136 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t155 + t192) - g(2) * (qJ(4) * t157 + t197) + t185 * t143 - t184 * qJD(4)) * MDP(11) + (-t106 * t181 - t129 * t130) * MDP(12) + (-t106 * t134 + t107 * t181 + t127 * t130 - t129 * t131) * MDP(13) + t108 * MDP(14) + t109 * MDP(15) + (qJD(3) * t127 + t137 * t107 + t120 * t134 + t124 * t131 - t182 * qJDD(5) + t189 * t154 + (qJD(4) * t181 - qJD(5) * t183) * qJD(5)) * MDP(17) + (qJD(3) * t129 + t137 * t106 - t120 * t181 - t124 * t130 - t183 * qJDD(5) + t189 * t156 + (qJD(4) * t134 + qJD(5) * t182) * qJD(5)) * MDP(18) + t217 * (t132 + t179 + t164); (-t114 * t167 + t115 * t169 - g(3)) * MDP(11) + t109 * MDP(17) - t108 * MDP(18) + (MDP(4) + MDP(7)) * (qJDD(2) - g(3)); (-qJD(1) * t138 + t187) * MDP(7) + (-qJD(1) * t136 - t178) * MDP(11) + (-qJD(1) * t127 + t108) * MDP(17) + (-qJD(1) * t129 + t109) * MDP(18) + (-MDP(6) - t217) * t175 + (MDP(7) * t149 + MDP(5) - t190) * qJDD(1); (qJD(1) * t184 + t132 + t189) * MDP(11) - t139 * MDP(17) + t140 * MDP(18) - t175 * t190 + ((t171 * MDP(17) + MDP(9)) * t169 + (t173 * MDP(17) - t171 * MDP(18) + MDP(8)) * t167) * qJDD(1) + ((t129 + t195) * MDP(17) + (-t127 + t180) * MDP(18)) * qJD(5); t129 * t127 * MDP(12) + (-t127 ^ 2 + t129 ^ 2) * MDP(13) + t186 * MDP(14) + (-t167 * t198 - t169 * t199 + t139) * MDP(15) + qJDD(5) * MDP(16) + (g(3) * t154 + t173 * t110 - t171 * t111 - t124 * t129 + t216 * t156) * MDP(17) + (g(3) * t156 - t171 * t110 - t173 * t111 + t124 * t127 - t216 * t154) * MDP(18) + ((t127 + t180) * MDP(14) + (t129 - t195) * MDP(15)) * qJD(5);];
tau = t1;
