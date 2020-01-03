% Calculate vector of inverse dynamics joint torques for
% S5RPRRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:26
% EndTime: 2020-01-03 11:52:28
% DurationCPUTime: 0.82s
% Computational Cost: add. (1092->150), mult. (1907->193), div. (0->0), fcn. (1095->14), ass. (0->93)
t189 = cos(pkin(9));
t175 = pkin(1) * t189 + pkin(2);
t166 = t175 * qJD(1);
t196 = cos(qJ(3));
t192 = sin(qJ(3));
t188 = sin(pkin(9));
t253 = pkin(1) * t188;
t229 = qJD(1) * t253;
t222 = t192 * t229;
t148 = t196 * t166 - t222;
t185 = qJD(1) + qJD(3);
t145 = pkin(3) * t185 + t148;
t191 = sin(qJ(4));
t149 = t166 * t192 + t196 * t229;
t195 = cos(qJ(4));
t239 = t149 * t195;
t136 = t145 * t191 + t239;
t184 = qJDD(1) + qJDD(3);
t164 = t175 * qJDD(1);
t235 = qJD(3) * t166;
t218 = t196 * t164 - t192 * t235;
t228 = qJD(1) * qJD(3) * t196;
t202 = (-qJDD(1) * t192 - t228) * t253 + t218;
t138 = pkin(3) * t184 + t202;
t141 = (qJDD(1) * t253 + t235) * t196 - qJD(3) * t222 + t164 * t192;
t183 = qJ(1) + pkin(9) + qJ(3);
t176 = qJ(4) + t183;
t170 = sin(t176);
t171 = cos(t176);
t209 = -g(2) * t171 - g(3) * t170 + t195 * t138 - t191 * t141;
t200 = -t136 * qJD(4) + t209;
t181 = qJDD(4) + t184;
t249 = pkin(4) * t181;
t215 = t249 + t200;
t234 = qJD(4) * t191;
t213 = g(2) * t170 - g(3) * t171 - t191 * t138 + t149 * t234;
t204 = t213 - (qJD(4) * t145 + t141) * t195;
t167 = t196 * t175;
t257 = -t192 * t253 + t167;
t154 = pkin(3) + t257;
t155 = t175 * t192 + t196 * t253;
t237 = t191 * t154 + t195 * t155;
t240 = t149 * t191;
t140 = t148 * t195 - t240;
t251 = pkin(3) * t191;
t177 = pkin(8) + t251;
t250 = pkin(3) * t195;
t178 = -pkin(4) - t250;
t182 = qJD(4) + t185;
t256 = -qJDD(5) * t177 + (-qJD(4) * t250 + t178 * t182 + t140) * qJD(5);
t173 = sin(t183);
t174 = cos(t183);
t255 = -g(2) * t174 - g(3) * t173;
t252 = pkin(3) * t182;
t152 = t257 * qJD(3);
t153 = t155 * qJD(3);
t243 = (t237 * qJD(4) + t152 * t191 + t153 * t195) * t182;
t242 = t136 * t182;
t241 = (t148 * t191 + t239) * t182;
t238 = qJDD(2) - g(1);
t190 = sin(qJ(5));
t186 = t190 ^ 2;
t194 = cos(qJ(5));
t236 = -t194 ^ 2 + t186;
t233 = qJD(5) * t182;
t232 = qJD(5) * t194;
t227 = -t145 - t252;
t135 = t145 * t195 - t240;
t133 = -pkin(4) * t182 - t135;
t221 = t133 * t232 - t215 * t190;
t193 = sin(qJ(1));
t197 = cos(qJ(1));
t219 = -g(2) * t197 - g(3) * t193;
t198 = qJD(5) ^ 2;
t162 = qJDD(5) * t190 + t194 * t198;
t163 = qJDD(5) * t194 - t190 * t198;
t217 = 0.2e1 * (t181 * t190 * t194 - t236 * t233) * MDP(12) + (0.2e1 * t182 * t190 * t232 + t181 * t186) * MDP(11) + t162 * MDP(13) + t163 * MDP(14) + t181 * MDP(8);
t216 = t154 * t195 - t155 * t191;
t212 = pkin(8) * t198 - t242 - t249;
t142 = -pkin(4) - t216;
t143 = pkin(8) + t237;
t211 = t142 * t181 + t143 * t198 + t243;
t210 = -pkin(8) * t181 - t133 * t182 + t204;
t129 = t216 * qJD(4) + t152 * t195 - t153 * t191;
t207 = -qJDD(5) * t143 + (t142 * t182 - t129) * qJD(5);
t206 = -pkin(4) * t233 - pkin(8) * qJDD(5) + qJD(5) * t135;
t205 = t177 * t198 + t178 * t181 + t234 * t252 - t241;
t131 = t133 * qJD(5) * t190;
t203 = t131 * MDP(16) + t221 * MDP(17) + t217;
t201 = g(2) * t173 - g(3) * t174 - t141;
t180 = t182 ^ 2;
t179 = t184 * MDP(5);
t1 = [qJDD(1) * MDP(1) + t219 * MDP(2) + (g(2) * t193 - g(3) * t197) * MDP(3) + (t219 + (t188 ^ 2 + t189 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t179 + (-t153 * t185 + t167 * t184 + (-t228 + (-qJDD(1) - t184) * t192) * t253 + t218 + t255) * MDP(6) + (-t152 * t185 - t155 * t184 + t201) * MDP(7) + (t216 * t181 + t200 - t243) * MDP(9) + (-t129 * t182 - t237 * t181 + t204) * MDP(10) + (t131 + t207 * t190 + (-t211 + t215) * t194) * MDP(16) + (t211 * t190 + t207 * t194 + t221) * MDP(17) + t217; t163 * MDP(16) - t162 * MDP(17) + t238 * MDP(4); t179 + (t149 * t185 + t202 + t255) * MDP(6) + (t148 * t185 + t201) * MDP(7) + (t181 * t250 + t209 + t241) * MDP(9) + (t140 * t182 - t195 * t141 - t181 * t251 + t213) * MDP(10) + (t227 * MDP(9) * t191 + (t227 * MDP(10) - t149 * MDP(9)) * t195) * qJD(4) + (t256 * MDP(16) + t205 * MDP(17)) * t190 + ((-t205 + t215) * MDP(16) + t256 * MDP(17)) * t194 + t203; (t200 + t242) * MDP(9) + (t135 * t182 + t204) * MDP(10) + (t206 * MDP(16) + t212 * MDP(17)) * t190 + ((-t212 + t215) * MDP(16) + t206 * MDP(17)) * t194 + t203; qJDD(5) * MDP(15) + t236 * MDP(12) * t180 + (t181 * MDP(14) + t238 * MDP(16) + t210 * MDP(17)) * t194 + (-t180 * t194 * MDP(11) + t181 * MDP(13) + t210 * MDP(16) - t238 * MDP(17)) * t190;];
tau = t1;
