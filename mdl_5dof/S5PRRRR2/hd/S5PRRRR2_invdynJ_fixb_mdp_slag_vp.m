% Calculate vector of inverse dynamics joint torques for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:55
% EndTime: 2019-12-05 17:04:57
% DurationCPUTime: 0.77s
% Computational Cost: add. (610->143), mult. (984->194), div. (0->0), fcn. (531->12), ass. (0->85)
t177 = cos(qJ(3));
t232 = pkin(2) * t177;
t158 = qJDD(2) * t232;
t166 = qJDD(2) + qJDD(3);
t173 = sin(qJ(3));
t230 = pkin(2) * qJD(2);
t208 = t173 * t230;
t131 = pkin(3) * t166 - qJD(3) * t208 + t158;
t167 = qJD(2) + qJD(3);
t139 = pkin(3) * t167 + t177 * t230;
t170 = qJ(2) + qJ(3);
t165 = qJ(4) + t170;
t155 = cos(t165);
t176 = cos(qJ(4));
t172 = sin(qJ(4));
t215 = qJD(4) * t172;
t238 = -g(2) * t155 + t176 * t131 - t139 * t215;
t213 = qJD(4) * t176;
t206 = t173 * t213;
t209 = qJDD(2) * t173;
t218 = qJD(3) * t177;
t237 = ((t172 * t218 + t206) * qJD(2) + t172 * t209) * pkin(2);
t236 = t237 - t238;
t187 = qJD(2) * t218 + t209;
t154 = sin(t165);
t201 = t172 * t208;
t190 = g(1) * t155 + g(2) * t154 + qJD(4) * t201 - t172 * t131;
t235 = -(t187 * pkin(2) + qJD(4) * t139) * t176 + t190;
t148 = g(1) * t154;
t128 = -t176 * t139 + t201;
t162 = qJD(4) + t167;
t229 = t128 * t162;
t228 = (t139 * t172 + t176 * t208) * t162;
t224 = t173 * t176;
t194 = t172 * t177 + t224;
t227 = t194 * t230 * t162;
t226 = t172 * MDP(9);
t225 = t172 * t173;
t223 = qJDD(1) - g(3);
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t212 = qJD(5) * t128;
t222 = t175 * t148 + t171 * t212;
t157 = pkin(3) + t232;
t221 = pkin(2) * t224 + t172 * t157;
t163 = sin(t170);
t164 = cos(t170);
t220 = g(1) * t164 + g(2) * t163;
t168 = t171 ^ 2;
t219 = -t175 ^ 2 + t168;
t216 = qJD(4) * t162;
t214 = qJD(4) * t173;
t211 = qJD(5) * t175;
t207 = t128 * t211 + t236 * t171;
t202 = qJD(3) * (-qJD(2) - t167);
t200 = g(1) * t163 - g(2) * t164 + t158;
t179 = qJD(5) ^ 2;
t198 = pkin(6) * t179 - t228;
t140 = qJDD(5) * t171 + t175 * t179;
t141 = qJDD(5) * t175 - t171 * t179;
t161 = qJDD(4) + t166;
t197 = 0.2e1 * (-t219 * t162 * qJD(5) + t161 * t171 * t175) * MDP(12) + (0.2e1 * t162 * t171 * t211 + t161 * t168) * MDP(11) + t140 * MDP(13) + t141 * MDP(14) + t161 * MDP(8);
t135 = pkin(2) * t225 - t157 * t176;
t196 = -(t157 * t215 + (t194 * qJD(3) + t206) * pkin(2)) * t162 - t135 * t161;
t156 = pkin(3) * t172 + pkin(6);
t195 = t156 * t179 - t227;
t193 = t176 * t177 - t225;
t192 = -pkin(6) * qJDD(5) - t212;
t134 = t193 * t230;
t191 = qJD(5) * t134 - qJDD(5) * t156;
t189 = -t176 * MDP(10) - MDP(7) - t226;
t188 = MDP(16) * t175 - MDP(17) * t171 + MDP(9);
t132 = pkin(6) + t221;
t185 = t132 * t179 - t196;
t184 = t148 + t238;
t183 = -pkin(6) * t161 - t229 + t235;
t121 = t157 * t213 + (t193 * qJD(3) - t172 * t214) * pkin(2);
t182 = -qJDD(5) * t132 + (t135 * t162 - t121) * qJD(5);
t181 = -t139 * t213 + t190;
t180 = t222 * MDP(16) + t207 * MDP(17) + t197;
t178 = cos(qJ(2));
t174 = sin(qJ(2));
t160 = t162 ^ 2;
t159 = t166 * MDP(5);
t1 = [t223 * MDP(1) + t141 * MDP(16) - t140 * MDP(17); qJDD(2) * MDP(2) + (g(1) * t174 - g(2) * t178) * MDP(3) + (g(1) * t178 + g(2) * t174) * MDP(4) + t159 + ((t166 * t177 + t173 * t202) * pkin(2) + t200) * MDP(6) + (((-qJDD(2) - t166) * t173 + t177 * t202) * pkin(2) + t220) * MDP(7) + (t184 + t196 - t237) * MDP(9) + (-t121 * t162 - t221 * t161 + t235) * MDP(10) + (t182 * t171 + (-t185 - t236) * t175 + t222) * MDP(16) + (t182 * t175 + (t185 - t148) * t171 + t207) * MDP(17) + t197; t159 + t200 * MDP(6) + t220 * MDP(7) + (t184 + t227) * MDP(9) + (t134 * t162 + t181) * MDP(10) + (t191 * MDP(16) + (t195 - t148) * MDP(17)) * t171 + ((-t195 - t236) * MDP(16) + t191 * MDP(17)) * t175 + (t189 * t209 + (((-qJD(3) + t167) * MDP(6) - MDP(9) * t213) * t173 + (t167 * MDP(7) + t189 * qJD(3)) * t177) * qJD(2)) * pkin(2) + ((-t161 * MDP(10) - t188 * t216) * t172 + (-MDP(10) * t216 + t188 * t161 + (MDP(16) * t171 + MDP(17) * t175) * qJD(5) * (-qJD(4) - t162)) * t176) * pkin(3) + t180; (t184 + t228) * MDP(9) + (t181 - t229) * MDP(10) + (t192 * MDP(16) + (t198 - t148) * MDP(17)) * t171 + (-t187 * t226 + (-qJD(2) * MDP(9) * t214 - t187 * MDP(10)) * t176) * pkin(2) + ((-t198 - t236) * MDP(16) + t192 * MDP(17)) * t175 + t180; qJDD(5) * MDP(15) + t219 * MDP(12) * t160 + (t161 * MDP(14) + t223 * MDP(16) + t183 * MDP(17)) * t175 + (-t160 * t175 * MDP(11) + t161 * MDP(13) + t183 * MDP(16) - t223 * MDP(17)) * t171;];
tau = t1;
