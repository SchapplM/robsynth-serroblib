% Calculate vector of inverse dynamics joint torques for
% S4RPRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:04
% EndTime: 2019-12-31 16:45:06
% DurationCPUTime: 1.06s
% Computational Cost: add. (733->179), mult. (1664->216), div. (0->0), fcn. (1135->8), ass. (0->86)
t237 = qJDD(1) * pkin(1);
t189 = sin(qJ(1));
t190 = cos(qJ(1));
t249 = g(1) * t189 - g(2) * t190;
t203 = -qJDD(2) + t237 + t249;
t210 = g(1) * t190 + g(2) * t189;
t185 = sin(pkin(6));
t186 = cos(pkin(6));
t248 = t186 * MDP(4) - t185 * MDP(5);
t188 = sin(qJ(3));
t242 = cos(qJ(3));
t160 = t185 * t242 + t188 * t186;
t155 = t160 * qJD(1);
t247 = qJ(2) * qJDD(1);
t221 = t242 * t186;
t231 = t188 * t185;
t200 = t221 - t231;
t238 = pkin(5) + qJ(2);
t165 = t238 * t185;
t166 = t238 * t186;
t201 = -t165 * t242 - t188 * t166;
t131 = qJD(2) * t200 + qJD(3) * t201;
t142 = -t188 * t165 + t166 * t242;
t184 = pkin(6) + qJ(3);
t179 = sin(t184);
t246 = -qJD(3) * t131 - qJDD(3) * t142 - t179 * t249;
t180 = cos(t184);
t226 = qJD(1) * qJD(2);
t244 = qJDD(1) * t238 + t226;
t147 = t244 * t185;
t148 = t244 * t186;
t161 = qJD(1) * t165;
t218 = qJD(3) * t242;
t223 = -t188 * t147 + t242 * t148 - t161 * t218;
t245 = -g(3) * t179 - t210 * t180 + t223;
t243 = t155 ^ 2;
t236 = qJDD(3) * pkin(3);
t211 = qJD(1) * t221;
t220 = qJD(1) * t231;
t153 = -t211 + t220;
t235 = t153 * t155;
t162 = qJD(1) * t166;
t232 = t188 * t162;
t230 = t185 ^ 2 + t186 ^ 2;
t140 = -t188 * t161 + t162 * t242;
t229 = qJD(3) * t140;
t228 = qJD(3) * t188;
t139 = -t161 * t242 - t232;
t227 = qJD(4) - t139;
t225 = qJDD(1) * t188;
t224 = qJDD(3) * qJ(4);
t216 = qJDD(1) * t242;
t222 = qJD(3) * t211 + t185 * t216 + t186 * t225;
t177 = pkin(2) * t186 + pkin(1);
t219 = t185 * t228;
t215 = t230 * qJD(1) ^ 2;
t214 = t139 + t232;
t213 = t242 * t147 + t188 * t148 - t161 * t228 + t162 * t218;
t212 = 0.2e1 * t230;
t208 = pkin(3) * t180 + qJ(4) * t179;
t206 = t185 * t225 - t186 * t216;
t202 = t177 + t208;
t164 = -qJD(1) * t177 + qJD(2);
t163 = -qJDD(1) * t177 + qJDD(2);
t158 = t160 * qJD(3);
t198 = -g(3) * t180 + t210 * t179 - t213;
t132 = qJD(2) * t160 + qJD(3) * t142;
t197 = -qJD(3) * t132 + qJDD(3) * t201 + t249 * t180;
t196 = t212 * t226 - t210;
t137 = qJD(1) * t219 - t222;
t138 = qJD(1) * t158 + t206;
t195 = pkin(3) * t138 + qJ(4) * t137 + t163;
t130 = pkin(3) * t153 - qJ(4) * t155 + t164;
t194 = t130 * t155 + qJDD(4) - t198;
t157 = -t186 * t218 + t219;
t152 = t153 ^ 2;
t136 = -pkin(3) * t200 - qJ(4) * t160 - t177;
t135 = pkin(3) * t155 + qJ(4) * t153;
t134 = qJD(3) * qJ(4) + t140;
t133 = -qJD(3) * pkin(3) + t227;
t129 = pkin(3) * t158 + qJ(4) * t157 - qJD(4) * t160;
t128 = (t153 - t220) * qJD(3) + t222;
t126 = qJDD(4) + t213 - t236;
t125 = t224 + (qJD(4) - t232) * qJD(3) + t223;
t124 = -qJD(4) * t155 + t195;
t1 = [qJDD(1) * MDP(1) + t249 * MDP(2) + t210 * MDP(3) + (t212 * t247 + t196) * MDP(6) + (pkin(1) * t203 + (t230 * t247 + t196) * qJ(2)) * MDP(7) + (-t137 * t160 - t155 * t157) * MDP(8) + (-t137 * t200 - t138 * t160 + t153 * t157 - t155 * t158) * MDP(9) + (-qJD(3) * t157 + qJDD(3) * t160) * MDP(10) + (-qJD(3) * t158 + qJDD(3) * t200) * MDP(11) + (-t138 * t177 + t158 * t164 - t163 * t200 + t197) * MDP(13) + (t137 * t177 - t157 * t164 + t160 * t163 + t246) * MDP(14) + (-t124 * t200 + t129 * t153 + t130 * t158 + t136 * t138 + t197) * MDP(15) + (t125 * t200 + t126 * t160 - t131 * t153 + t132 * t155 - t133 * t157 - t134 * t158 + t137 * t201 - t138 * t142 - t210) * MDP(16) + (-t124 * t160 - t129 * t155 + t130 * t157 + t136 * t137 - t246) * MDP(17) + (t124 * t136 + t125 * t142 - t126 * t201 + t130 * t129 + t134 * t131 + t133 * t132 + (-g(1) * t238 - g(2) * t202) * t190 + (g(1) * t202 - g(2) * t238) * t189) * MDP(18) + t248 * (t203 + t237); -MDP(6) * t215 + (-qJ(2) * t215 - t203) * MDP(7) + (-t152 - t243) * MDP(16) + (t134 * t153 + (-qJD(4) - t133) * t155 + t195 - t249) * MDP(18) + (MDP(13) + MDP(15)) * (0.2e1 * t155 * qJD(3) + t206) + (-MDP(14) + MDP(17)) * ((t153 + t220) * qJD(3) - t222) - t248 * qJDD(1); MDP(8) * t235 + (-t152 + t243) * MDP(9) + t128 * MDP(10) - t206 * MDP(11) + qJDD(3) * MDP(12) + (-t155 * t164 + t198 + t229) * MDP(13) + (qJD(3) * t214 + t153 * t164 - t245) * MDP(14) + (-t135 * t153 - t194 + t229 + 0.2e1 * t236) * MDP(15) + (pkin(3) * t137 - qJ(4) * t138 + (t134 - t140) * t155 + (t133 - t227) * t153) * MDP(16) + (0.2e1 * t224 - t130 * t153 + t135 * t155 + (0.2e1 * qJD(4) - t214) * qJD(3) + t245) * MDP(17) + (-t126 * pkin(3) - g(3) * t208 + t125 * qJ(4) - t130 * t135 - t133 * t140 + t134 * t227 + t210 * (pkin(3) * t179 - qJ(4) * t180)) * MDP(18); (-qJDD(3) + t235) * MDP(15) + t128 * MDP(16) + (-qJD(3) ^ 2 - t243) * MDP(17) + (-qJD(3) * t134 + t194 - t236) * MDP(18);];
tau = t1;
