% Calculate vector of inverse dynamics joint torques for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:32
% EndTime: 2021-01-15 14:56:34
% DurationCPUTime: 0.82s
% Computational Cost: add. (512->183), mult. (1041->223), div. (0->0), fcn. (636->6), ass. (0->94)
t161 = cos(qJ(4));
t160 = sin(qJ(3));
t140 = qJD(3) * pkin(6) + t160 * qJD(2);
t182 = qJ(5) * qJD(3) + t140;
t226 = t182 * t161;
t162 = cos(qJ(3));
t194 = qJDD(2) * t160;
t196 = qJD(2) * qJD(3);
t129 = qJDD(3) * pkin(6) + t162 * t196 + t194;
t179 = -qJ(5) * qJDD(3) - t129;
t169 = qJD(3) * qJD(5) - t179;
t159 = sin(qJ(4));
t195 = qJD(3) * qJD(4);
t184 = t159 * t195;
t197 = qJD(1) * qJD(4);
t198 = qJD(4) * t159;
t186 = -t159 * qJDD(1) - t140 * t198 - t161 * t197;
t114 = -qJ(5) * t184 + t169 * t161 + t186;
t119 = -qJD(1) * t161 - t182 * t159;
t217 = qJD(4) * pkin(4);
t118 = t119 + t217;
t225 = -qJD(4) * t118 + t114;
t155 = t159 ^ 2;
t156 = t161 ^ 2;
t206 = t155 + t156;
t224 = t206 * MDP(15);
t223 = pkin(4) * t155;
t215 = sin(pkin(7));
t216 = cos(pkin(7));
t130 = -t215 * t160 - t216 * t162;
t222 = g(1) * t130;
t131 = t216 * t160 - t215 * t162;
t221 = g(1) * t131;
t220 = g(2) * t130;
t219 = g(2) * t131;
t154 = g(3) * t161;
t158 = qJ(5) + pkin(6);
t218 = qJD(3) * pkin(3);
t213 = qJDD(4) * pkin(4);
t212 = t131 * t159;
t164 = qJD(3) ^ 2;
t211 = t161 * t164;
t210 = t118 - t119;
t146 = t160 * t196;
t204 = qJD(2) * t162;
t185 = qJD(4) * t204;
t209 = t161 * t146 + t159 * t185;
t208 = -g(1) * t212 + t161 * t185;
t207 = t155 - t156;
t163 = qJD(4) ^ 2;
t205 = t163 + t164;
t149 = pkin(4) * t161 + pkin(3);
t202 = qJD(3) * t149;
t127 = qJD(5) - t202 - t204;
t203 = qJD(3) * t127;
t201 = qJD(3) * t159;
t200 = qJD(3) * t161;
t193 = qJDD(2) * t162;
t192 = qJDD(3) * t149;
t191 = qJDD(3) * t159;
t190 = qJDD(3) * t161;
t189 = qJDD(4) * t159;
t188 = MDP(11) + MDP(13);
t187 = MDP(12) + MDP(14);
t183 = qJD(4) * t158;
t141 = -t204 - t218;
t181 = -qJD(3) * t141 - t129;
t180 = 0.2e1 * t161 * t195;
t178 = t146 - t193;
t142 = pkin(4) * t184;
t177 = -t161 * t222 - t186;
t176 = -t220 + t221;
t175 = t219 + t222;
t174 = -g(1) * t215 + g(2) * t216;
t151 = t159 * qJD(1);
t120 = -t151 + t226;
t173 = t118 * t159 - t120 * t161;
t172 = -MDP(5) + t224;
t117 = qJDD(5) + t142 + t178 - t192;
t171 = t117 - t192 + t220;
t170 = 0.2e1 * t184 - t190;
t147 = t159 * t197;
t168 = -g(2) * t212 - t161 * qJDD(1) - t159 * t222 + t147 + t154;
t167 = 0.2e1 * qJDD(3) * pkin(3) - pkin(6) * t163 - t178 - t220;
t166 = -pkin(6) * qJDD(4) + (t141 - t218) * qJD(4);
t165 = (-qJD(5) - t127) * qJD(3) + t179;
t136 = t158 * t161;
t135 = t158 * t159;
t134 = qJDD(4) * t161 - t163 * t159;
t133 = t161 * t163 + t189;
t126 = -qJD(5) * t159 - t161 * t183;
t125 = qJD(5) * t161 - t159 * t183;
t113 = t213 + t147 + (-t182 * qJD(4) - qJDD(1)) * t161 - t169 * t159;
t1 = [(t173 * qJD(4) - t113 * t161 - t114 * t159 - g(3)) * MDP(16) + (MDP(1) + MDP(2)) * (qJDD(1) - g(3)) - t188 * t134 + t187 * t133; (qJDD(2) + t174) * MDP(2) + t174 * MDP(16) + t188 * (-t170 * t162 + (-t205 * t161 - t189) * t160) + t187 * ((-qJDD(4) * t160 - 0.2e1 * t162 * t195) * t161 + (-qJDD(3) * t162 + t205 * t160) * t159) + (qJDD(3) * MDP(4) + (-t118 * t201 + t120 * t200 - t117) * MDP(16) + t172 * t164) * t162 + (-t164 * MDP(4) + (-t113 * t159 - t120 * t198 + t225 * t161 + t203) * MDP(16) + t172 * qJDD(3)) * t160; qJDD(3) * MDP(3) + (t176 + t193) * MDP(4) + (-t175 - t194) * MDP(5) + (qJDD(3) * t155 + t159 * t180) * MDP(6) + 0.2e1 * (t159 * t190 - t207 * t195) * MDP(7) + t133 * MDP(8) + t134 * MDP(9) + (t166 * t159 + (t167 + t221) * t161 + t209) * MDP(11) + (t166 * t161 + (-t167 - t146) * t159 + t208) * MDP(12) + (-qJDD(4) * t135 + (t126 + (t127 - t202) * t159) * qJD(4) + (-t142 - t171 + t221) * t161 + t209) * MDP(13) + (-qJDD(4) * t136 + (t171 - t146) * t159 + (t127 * t161 - t125 + (-t149 * t161 + t223) * qJD(3)) * qJD(4) + t208) * MDP(14) + ((qJDD(3) * t136 + t225) * t161 + (-qJD(4) * t120 + qJDD(3) * t135 - t113) * t159 + (t125 * t161 - t126 * t159 + (t135 * t161 - t136 * t159) * qJD(4) - t206 * t204) * qJD(3) + t175) * MDP(15) + (t114 * t136 + t120 * t125 - t113 * t135 + t118 * t126 - t117 * t149 + t127 * pkin(4) * t198 - g(1) * (-t130 * t158 - t131 * t149) - g(2) * (t130 * t149 - t131 * t158) + (-t127 * t160 + t173 * t162) * qJD(2)) * MDP(16); -t159 * MDP(6) * t211 + t207 * t164 * MDP(7) + MDP(8) * t191 + MDP(9) * t190 + qJDD(4) * MDP(10) + (-t151 * qJD(4) + t181 * t159 + t168) * MDP(11) + ((-qJD(4) * t140 - g(3)) * t159 + (t181 - t197 - t219) * t161 + t177) * MDP(12) + (0.2e1 * t213 + (t120 - t226) * qJD(4) + (pkin(4) * t211 + t165) * t159 + t168) * MDP(13) + (-t164 * t223 - g(3) * t159 + (qJ(5) * t201 + t119) * qJD(4) + (t165 - t219) * t161 + t177) * MDP(14) + (-pkin(4) * t191 + (t210 - t217) * t200) * MDP(15) + (t210 * t120 + (t154 + t113 + (-t175 - t203) * t159) * pkin(4)) * MDP(16); t170 * MDP(13) + (t180 + t191) * MDP(14) + (t173 * qJD(3) + t117 - t176) * MDP(16) - t164 * t224;];
tau = t1;
