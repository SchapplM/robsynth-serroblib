% Calculate vector of inverse dynamics joint torques for
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:17
% EndTime: 2019-12-31 17:22:18
% DurationCPUTime: 0.75s
% Computational Cost: add. (633->149), mult. (951->199), div. (0->0), fcn. (505->12), ass. (0->85)
t177 = cos(qJ(2));
t231 = pkin(1) * t177;
t158 = qJDD(1) * t231;
t166 = qJDD(1) + qJDD(2);
t173 = sin(qJ(2));
t228 = pkin(1) * qJD(1);
t207 = t173 * t228;
t133 = pkin(2) * t166 - qJD(2) * t207 + t158;
t167 = qJD(1) + qJD(2);
t140 = pkin(2) * t167 + t177 * t228;
t170 = qJ(1) + qJ(2);
t165 = qJ(3) + t170;
t154 = cos(t165);
t176 = cos(qJ(3));
t172 = sin(qJ(3));
t214 = qJD(3) * t172;
t236 = -g(2) * t154 + t176 * t133 - t140 * t214;
t212 = qJD(3) * t176;
t204 = t173 * t212;
t209 = qJDD(1) * t173;
t216 = qJD(2) * t177;
t161 = qJDD(3) + t166;
t230 = pkin(3) * t161;
t235 = -t230 + (t172 * t209 + (t172 * t216 + t204) * qJD(1)) * pkin(1) - t236;
t203 = qJD(1) * t216;
t190 = t203 + t209;
t153 = sin(t165);
t198 = t172 * t207;
t192 = g(1) * t154 + g(2) * t153 + qJD(3) * t198 - t172 * t133;
t234 = -(t190 * pkin(1) + qJD(3) * t140) * t176 + t192;
t222 = t172 * t173;
t193 = t176 * t177 - t222;
t137 = t193 * t228;
t155 = pkin(2) * t172 + pkin(7);
t156 = -pkin(2) * t176 - pkin(3);
t162 = qJD(3) + t167;
t233 = -qJDD(4) * t155 + (-pkin(2) * t212 + t156 * t162 + t137) * qJD(4);
t147 = g(1) * t153;
t227 = MDP(8) * t172;
t157 = pkin(2) + t231;
t221 = t173 * t176;
t194 = t172 * t177 + t221;
t226 = (t157 * t214 + (t194 * qJD(2) + t204) * pkin(1)) * t162;
t225 = (t140 * t172 + t176 * t207) * t162;
t224 = t194 * t228 * t162;
t223 = t161 * t176;
t130 = t140 * t176 - t198;
t128 = -pkin(3) * t162 - t130;
t171 = sin(qJ(4));
t175 = cos(qJ(4));
t220 = t128 * qJD(4) * t171 + t175 * t147;
t219 = pkin(1) * t221 + t172 * t157;
t163 = sin(t170);
t164 = cos(t170);
t218 = g(1) * t164 + g(2) * t163;
t168 = t171 ^ 2;
t217 = -t175 ^ 2 + t168;
t213 = qJD(3) * t173;
t211 = qJD(4) * t162;
t210 = qJD(4) * t175;
t206 = t128 * t210 + t235 * t171;
t205 = t162 * t214;
t202 = qJD(1) * t213;
t199 = qJD(2) * (-qJD(1) - t167);
t197 = g(1) * t163 - g(2) * t164 + t158;
t179 = qJD(4) ^ 2;
t195 = 0.2e1 * (t171 * t161 * t175 - t217 * t211) * MDP(11) + (0.2e1 * t162 * t171 * t210 + t161 * t168) * MDP(10) + (qJDD(4) * t171 + t175 * t179) * MDP(12) + (qJDD(4) * t175 - t171 * t179) * MDP(13) + t161 * MDP(7);
t191 = -MDP(9) * t176 - MDP(6) - t227;
t189 = pkin(7) * t179 - t225 - t230;
t134 = pkin(1) * t222 - t157 * t176 - pkin(3);
t135 = pkin(7) + t219;
t188 = t134 * t161 + t135 * t179 + t226;
t187 = t147 + t236;
t186 = -pkin(7) * t161 - t128 * t162 + t234;
t123 = t157 * t212 + (qJD(2) * t193 - t172 * t213) * pkin(1);
t184 = -qJDD(4) * t135 + (t134 * t162 - t123) * qJD(4);
t183 = -pkin(3) * t211 - pkin(7) * qJDD(4) + qJD(4) * t130;
t182 = -t140 * t212 + t192;
t181 = pkin(2) * t205 + t155 * t179 + t156 * t161 - t224;
t180 = MDP(15) * t220 + MDP(16) * t206 + t195;
t178 = cos(qJ(1));
t174 = sin(qJ(1));
t160 = t162 ^ 2;
t159 = t166 * MDP(4);
t1 = [qJDD(1) * MDP(1) + (g(1) * t174 - g(2) * t178) * MDP(2) + (g(1) * t178 + g(2) * t174) * MDP(3) + t159 + ((t166 * t177 + t173 * t199) * pkin(1) + t197) * MDP(5) + (((-qJDD(1) - t166) * t173 + t177 * t199) * pkin(1) + t218) * MDP(6) + (t157 * t223 - t226 + (-t176 * t202 + (-t203 + (-qJDD(1) - t161) * t173) * t172) * pkin(1) + t187) * MDP(8) + (-t123 * t162 - t219 * t161 + t234) * MDP(9) + (t184 * t171 + (-t188 - t235) * t175 + t220) * MDP(15) + (t184 * t175 + (t188 - t147) * t171 + t206) * MDP(16) + t195; t159 + t197 * MDP(5) + t218 * MDP(6) + (t187 + t224) * MDP(8) + (t137 * t162 + t182) * MDP(9) + ((-t205 + t223) * MDP(8) + (-t161 * t172 - t162 * t212) * MDP(9)) * pkin(2) + (t191 * t209 + (((-qJD(2) + t167) * MDP(5) - MDP(8) * t212) * t173 + (t167 * MDP(6) + qJD(2) * t191) * t177) * qJD(1)) * pkin(1) + ((t181 - t147) * MDP(16) + t233 * MDP(15)) * t171 + ((-t181 - t235) * MDP(15) + t233 * MDP(16)) * t175 + t180; (t187 + t225) * MDP(8) + (t130 * t162 + t182) * MDP(9) + (-t190 * t227 + (-MDP(8) * t202 - MDP(9) * t190) * t176) * pkin(1) + (t183 * MDP(15) + (t189 - t147) * MDP(16)) * t171 + ((-t189 - t235) * MDP(15) + t183 * MDP(16)) * t175 + t180; qJDD(4) * MDP(14) + t217 * MDP(11) * t160 + (t161 * MDP(13) - g(3) * MDP(15) + MDP(16) * t186) * t175 + (-t160 * t175 * MDP(10) + t161 * MDP(12) + MDP(15) * t186 + g(3) * MDP(16)) * t171;];
tau = t1;
