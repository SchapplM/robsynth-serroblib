% Calculate vector of inverse dynamics joint torques for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:52
% EndTime: 2019-12-31 18:17:54
% DurationCPUTime: 0.72s
% Computational Cost: add. (620->139), mult. (960->168), div. (0->0), fcn. (502->12), ass. (0->81)
t166 = cos(pkin(8));
t153 = pkin(1) * t166 + pkin(2);
t168 = sin(qJ(3));
t165 = sin(pkin(8));
t218 = pkin(1) * t165;
t197 = qJD(3) * t218;
t171 = cos(qJ(3));
t206 = qJD(3) * t171;
t186 = t153 * t206 - t168 * t197;
t122 = -qJD(4) - t186;
t185 = t153 * t168 + t171 * t218;
t126 = qJ(4) + t185;
t159 = qJDD(1) + qJDD(3);
t160 = qJD(1) + qJD(3);
t161 = qJ(1) + pkin(8);
t157 = qJ(3) + t161;
t151 = sin(t157);
t152 = cos(t157);
t221 = g(1) * t152 + g(2) * t151;
t223 = -t122 * t160 + t126 * t159 - t221;
t222 = -qJD(1) * t197 + t153 * qJDD(1);
t220 = -g(1) * t151 + g(2) * t152;
t198 = qJD(1) * t218;
t140 = t153 * qJD(1);
t212 = t140 * t168;
t121 = t171 * t198 + t212;
t214 = qJ(4) * t160;
t119 = t121 + t214;
t219 = -t119 * t160 + t220;
t158 = t160 ^ 2;
t173 = -pkin(3) - pkin(7);
t217 = pkin(3) * t159;
t169 = sin(qJ(1));
t216 = t169 * pkin(1);
t172 = cos(qJ(1));
t215 = t172 * pkin(1);
t154 = t159 * qJ(4);
t211 = qJDD(2) - g(3);
t156 = t160 * qJD(4);
t196 = qJDD(1) * t218;
t182 = -t140 * t206 - t222 * t168 - t171 * t196;
t113 = -t154 - t156 + t182;
t167 = sin(qJ(5));
t170 = cos(qJ(5));
t204 = qJD(5) * t170;
t210 = -t113 * t167 + t119 * t204;
t209 = t152 * pkin(3) + t151 * qJ(4);
t174 = qJD(5) ^ 2;
t208 = -t158 - t174;
t163 = t170 ^ 2;
t207 = t167 ^ 2 - t163;
t205 = qJD(5) * t160;
t120 = -t171 * t140 + t168 * t198;
t203 = qJD(4) + t120;
t193 = t153 * t171 - t168 * t218;
t127 = -pkin(3) - t193;
t125 = -pkin(7) + t127;
t202 = qJDD(5) * t125;
t201 = qJDD(5) * t167;
t200 = qJDD(5) * t170;
t199 = qJDD(5) * t173;
t195 = -t151 * pkin(3) + t152 * qJ(4);
t194 = -t121 + t214;
t191 = -qJD(3) * t212 - t168 * t196 + t222 * t171;
t136 = -t170 * t174 - t201;
t137 = -t167 * t174 + t200;
t189 = 0.2e1 * (-t159 * t167 * t170 + t207 * t205) * MDP(12) + (-0.2e1 * t160 * t167 * t204 + t159 * t163) * MDP(11) + t136 * MDP(14) + t137 * MDP(13) + t159 * MDP(5);
t187 = qJDD(4) - t191;
t184 = -t191 + t220;
t114 = t187 - t217;
t183 = -t173 * t159 - t187 - t219;
t181 = t121 * t160 - t184;
t124 = t185 * qJD(3);
t180 = -t124 * t160 - t184;
t179 = t182 + t221;
t178 = -t125 * t174 + t223;
t177 = -t120 * t160 + t179;
t176 = t203 * t160 - t173 * t174 + t154 - t221;
t118 = -pkin(3) * t160 + t203;
t111 = t113 * t170;
t1 = [(g(1) * t169 - g(2) * t172) * MDP(2) + (g(1) * t172 + g(2) * t169) * MDP(3) + (g(1) * t216 - g(2) * t215) * MDP(4) + (t193 * t159 + t180) * MDP(6) + (-t185 * t159 - t186 * t160 + t179) * MDP(7) + (qJDD(4) + (-pkin(3) + t127) * t159 - t180) * MDP(8) + (-t113 + t223) * MDP(9) + (-t113 * t126 - t119 * t122 + t114 * t127 + t118 * t124 - g(1) * (-pkin(2) * sin(t161) - t216 + t195) - g(2) * (pkin(2) * cos(t161) + t215 + t209)) * MDP(10) + t210 * MDP(16) - t111 * MDP(17) + (MDP(1) + (t165 ^ 2 + t166 ^ 2) * MDP(4) * pkin(1) ^ 2) * qJDD(1) + ((qJD(5) * t124 + t126 * t205 + t202) * MDP(16) + t178 * MDP(17)) * t170 + (t178 * MDP(16) + (-t202 + (-t126 * t160 - t119 - t124) * qJD(5)) * MDP(17)) * t167 + t189; MDP(16) * t136 - MDP(17) * t137 + (MDP(4) + MDP(10)) * t211; t181 * MDP(6) + t177 * MDP(7) + (qJDD(4) - t181 - 0.2e1 * t217) * MDP(8) + (0.2e1 * t154 + 0.2e1 * t156 - t177) * MDP(9) + (-t114 * pkin(3) - g(1) * t195 - g(2) * t209 - t113 * qJ(4) - t118 * t121 + t203 * t119) * MDP(10) + ((t194 * qJD(5) + t199) * t170 + t176 * t167 + t210) * MDP(16) + (-t111 + (-t199 + (-t119 - t194) * qJD(5)) * t167 + t176 * t170) * MDP(17) + t189; t159 * MDP(8) - t158 * MDP(9) + (t114 + t219) * MDP(10) + (t208 * t167 + t200) * MDP(16) + (t208 * t170 - t201) * MDP(17); qJDD(5) * MDP(15) - t207 * MDP(12) * t158 + (t159 * MDP(13) - t183 * MDP(16) - t211 * MDP(17)) * t170 + (t170 * t158 * MDP(11) - t159 * MDP(14) - t211 * MDP(16) + t183 * MDP(17)) * t167;];
tau = t1;
