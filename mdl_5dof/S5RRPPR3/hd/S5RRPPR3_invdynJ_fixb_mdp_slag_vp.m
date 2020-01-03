% Calculate vector of inverse dynamics joint torques for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:38
% EndTime: 2019-12-31 19:26:40
% DurationCPUTime: 0.77s
% Computational Cost: add. (658->169), mult. (993->205), div. (0->0), fcn. (531->12), ass. (0->98)
t177 = qJD(1) + qJD(2);
t188 = cos(qJ(2));
t233 = pkin(1) * qJD(1);
t209 = t188 * t233;
t145 = pkin(2) * t177 + t209;
t182 = sin(pkin(8));
t183 = cos(pkin(8));
t185 = sin(qJ(2));
t210 = t185 * t233;
t132 = t145 * t182 + t183 * t210;
t130 = qJ(4) * t177 + t132;
t181 = qJ(1) + qJ(2);
t171 = pkin(8) + t181;
t161 = sin(t171);
t162 = cos(t171);
t238 = pkin(1) * t188;
t169 = qJDD(1) * t238;
t176 = qJDD(1) + qJDD(2);
t232 = pkin(1) * qJD(2);
t208 = qJD(1) * t232;
t137 = pkin(2) * t176 - t185 * t208 + t169;
t215 = qJDD(1) * t185;
t241 = pkin(1) * t215 + t188 * t208;
t124 = t137 * t183 - t182 * t241;
t198 = qJDD(4) - t124;
t193 = -g(1) * t161 + g(2) * t162 + t198;
t242 = -t130 * t177 + t193;
t172 = sin(t181);
t165 = g(1) * t172;
t173 = cos(t181);
t240 = -g(2) * t173 + t165;
t175 = t177 ^ 2;
t239 = pkin(1) * t185;
t237 = pkin(3) * t176;
t186 = sin(qJ(1));
t235 = g(1) * t186;
t231 = MDP(7) * t183;
t158 = t182 * t239;
t153 = qJD(2) * t158;
t135 = t183 * t188 * t232 + qJD(4) - t153;
t229 = t135 * t177;
t168 = pkin(2) + t238;
t228 = t168 * t182;
t227 = t183 * t185;
t226 = qJDD(3) - g(3);
t125 = t182 * t137 + t183 * t241;
t221 = qJD(4) * t177;
t122 = qJ(4) * t176 + t125 + t221;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t219 = qJD(5) * t187;
t225 = t122 * t184 + t130 * t219;
t224 = g(1) * t173 + g(2) * t172;
t190 = qJD(5) ^ 2;
t223 = -t175 - t190;
t179 = t187 ^ 2;
t222 = t184 ^ 2 - t179;
t220 = qJD(5) * t177;
t218 = MDP(7) + MDP(10);
t217 = -qJD(1) - t177;
t152 = t182 * t210;
t142 = t183 * t209 - t152;
t216 = qJD(4) - t142;
t204 = t168 * t183 - t158;
t139 = -pkin(3) - t204;
t136 = -pkin(7) + t139;
t214 = qJDD(5) * t136;
t163 = -pkin(2) * t183 - pkin(3);
t159 = -pkin(7) + t163;
t213 = qJDD(5) * t159;
t212 = qJDD(5) * t184;
t211 = qJDD(5) * t187;
t167 = pkin(2) * t173;
t206 = t162 * pkin(3) + t161 * qJ(4) + t167;
t197 = pkin(1) * (t182 * t188 + t227);
t140 = qJD(1) * t197;
t160 = pkin(2) * t182 + qJ(4);
t205 = t160 * t177 - t140;
t131 = t145 * t183 - t152;
t203 = qJD(1) * (-qJD(2) + t177);
t201 = t169 + t240;
t200 = -g(1) * t162 - g(2) * t161;
t146 = -t187 * t190 - t212;
t147 = -t184 * t190 + t211;
t199 = 0.2e1 * (-t176 * t184 * t187 + t220 * t222) * MDP(12) + (-0.2e1 * t177 * t184 * t219 + t176 * t179) * MDP(11) + t146 * MDP(14) + t147 * MDP(13) + t176 * MDP(4);
t196 = g(1) * (-pkin(2) * t172 - pkin(3) * t161 + t162 * qJ(4));
t195 = t200 + t125;
t194 = -(-pkin(3) - pkin(7)) * t176 - t242;
t138 = pkin(1) * t227 + qJ(4) + t228;
t192 = -t136 * t190 + t138 * t176 + t200 + t229;
t191 = -t159 * t190 + t160 * t176 + t177 * t216 + t200;
t189 = cos(qJ(1));
t174 = t189 * pkin(1);
t141 = qJD(2) * t197;
t129 = -pkin(3) * t177 + qJD(4) - t131;
t123 = t198 - t237;
t120 = t122 * t187;
t1 = [((-pkin(3) + t139) * MDP(8) + (qJ(4) + t138) * MDP(9)) * t176 + t225 * MDP(16) + t120 * MDP(17) + (t125 * t228 - t132 * t153 + t124 * t204 - t131 * t141 + pkin(2) * t165 - g(2) * (t167 + t174)) * MDP(7) + ((qJD(5) * t141 + t138 * t220 + t214) * MDP(16) + t192 * MDP(17)) * t187 + (t192 * MDP(16) + (-t214 + (-t138 * t177 - t130 - t141) * qJD(5)) * MDP(17)) * t184 + (t141 * t177 + t193) * MDP(8) + (t195 + t221 + t229) * MDP(9) + (t122 * t138 + t130 * t135 + t123 * t139 + t129 * t141 - t196 - g(2) * (t174 + t206)) * MDP(10) + qJDD(1) * MDP(1) + (t176 * t188 * MDP(5) + t218 * t235 + ((-qJDD(1) - t176) * MDP(6) + t125 * t231) * t185 + (t217 * MDP(5) * t185 + (MDP(6) * t217 + t132 * t231) * t188) * qJD(2)) * pkin(1) + t201 * MDP(5) + t224 * MDP(6) + t199 + (-g(2) * t189 + t235) * MDP(2) + (g(1) * t189 + g(2) * t186) * MDP(3); (t203 * t239 + t201) * MDP(5) + ((t188 * t203 - t215) * pkin(1) + t224) * MDP(6) + (t131 * t140 - t132 * t142 + (t124 * t183 + t125 * t182 + t240) * pkin(2)) * MDP(7) + (-t140 * t177 + (-pkin(3) + t163) * t176 + t193) * MDP(8) + ((0.2e1 * qJD(4) - t142) * t177 + (qJ(4) + t160) * t176 + t195) * MDP(9) + (-g(2) * t206 + t122 * t160 + t123 * t163 - t129 * t140 + t130 * t216 - t196) * MDP(10) + ((qJD(5) * t205 + t213) * t187 + t191 * t184 + t225) * MDP(16) + (t120 + (-t213 + (-t130 - t205) * qJD(5)) * t184 + t191 * t187) * MDP(17) + t199; MDP(16) * t146 - MDP(17) * t147 + t218 * t226; t176 * MDP(8) - t175 * MDP(9) + (-t237 + t242) * MDP(10) + (t184 * t223 + t211) * MDP(16) + (t187 * t223 - t212) * MDP(17); qJDD(5) * MDP(15) - t222 * MDP(12) * t175 + (t176 * MDP(13) - MDP(16) * t194 - MDP(17) * t226) * t187 + (t187 * t175 * MDP(11) - t176 * MDP(14) - MDP(16) * t226 + MDP(17) * t194) * t184;];
tau = t1;
