% Calculate vector of inverse dynamics joint torques for
% S4RPRR9
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:32
% EndTime: 2019-12-31 16:56:34
% DurationCPUTime: 1.17s
% Computational Cost: add. (555->191), mult. (1097->272), div. (0->0), fcn. (653->6), ass. (0->93)
t149 = sin(qJ(1));
t152 = cos(qJ(1));
t169 = g(1) * t149 - g(2) * t152;
t222 = qJDD(2) - t169;
t151 = cos(qJ(3));
t181 = qJDD(1) * t151;
t221 = qJD(3) * qJD(4) + t181;
t214 = pkin(1) * qJDD(1);
t220 = t214 - t222;
t153 = -pkin(1) - pkin(5);
t135 = t153 * qJDD(1) + qJDD(2);
t136 = t153 * qJD(1) + qJD(2);
t148 = sin(qJ(3));
t191 = qJD(3) * t148;
t118 = -qJDD(3) * pkin(3) - t135 * t151 + t136 * t191;
t207 = t136 * t151;
t126 = -qJD(3) * pkin(3) - t207;
t185 = qJD(1) * qJD(3);
t177 = t151 * t185;
t182 = qJDD(1) * t148;
t127 = qJDD(4) + t177 + t182;
t133 = pkin(3) * t148 - pkin(6) * t151 + qJ(2);
t138 = qJD(1) * t148 + qJD(4);
t120 = t133 * qJD(1);
t190 = qJD(3) * t151;
t172 = -qJDD(3) * pkin(6) - qJD(4) * t120 - t135 * t148 - t136 * t190;
t219 = -(qJD(4) * t133 + t153 * t190) * t138 + t118 * t151 + (-qJD(3) * t126 - t127 * t153 + t172) * t148;
t171 = pkin(3) * t151 + pkin(6) * t148;
t216 = g(3) * t148;
t218 = (pkin(6) * qJD(4) + t171 * qJD(1)) * t138 + t169 * t151 + t118 - t216;
t215 = g(3) * t151;
t155 = qJD(1) ^ 2;
t213 = qJ(2) * t155;
t147 = sin(qJ(4));
t150 = cos(qJ(4));
t186 = t150 * qJD(3);
t187 = qJD(4) * t151;
t162 = -t147 * t187 - t148 * t186;
t180 = t147 * qJDD(3) + t221 * t150;
t114 = t162 * qJD(1) + t180;
t212 = t114 * t147;
t195 = qJD(1) * t151;
t129 = t147 * t195 - t186;
t211 = t129 * t138;
t192 = qJD(3) * t147;
t131 = t150 * t195 + t192;
t210 = t131 * t138;
t209 = t133 * t127;
t208 = t136 * t148;
t206 = t147 * t148;
t205 = t147 * t149;
t204 = t147 * t152;
t203 = t148 * t153;
t202 = t149 * t150;
t201 = t150 * t138;
t200 = t150 * t152;
t199 = t151 * t153;
t198 = -t150 * qJDD(3) - t185 * t206;
t146 = t151 ^ 2;
t197 = t148 ^ 2 - t146;
t154 = qJD(3) ^ 2;
t196 = -t154 - t155;
t194 = qJD(3) * t129;
t193 = qJD(3) * t131;
t125 = qJD(3) * pkin(6) + t208;
t189 = qJD(4) * t125;
t188 = qJD(4) * t138;
t183 = qJDD(1) * qJ(2);
t179 = t138 * t192;
t178 = t138 * t186;
t128 = t171 * qJD(3) + qJD(2);
t117 = t128 * qJD(1) + t133 * qJDD(1);
t173 = -t117 + t189;
t170 = g(1) * t152 + g(2) * t149;
t168 = t172 + t215;
t166 = t147 * t127 + t150 * t188;
t165 = t150 * t127 - t147 * t188;
t163 = 0.2e1 * qJ(2) * t185 + qJDD(3) * t153;
t161 = 0.2e1 * qJD(1) * qJD(2) - t170;
t160 = -t135 + t169 + t213;
t159 = t161 + 0.2e1 * t183;
t158 = -pkin(6) * t127 + (t126 + t207) * t138;
t157 = -t153 * t154 + t159;
t142 = qJDD(3) * t151;
t124 = t148 * t200 - t205;
t123 = t148 * t204 + t202;
t122 = t148 * t202 + t204;
t121 = -t148 * t205 + t200;
t116 = t150 * t117;
t115 = t131 * qJD(4) + t147 * t181 + t198;
t113 = t120 * t147 + t125 * t150;
t112 = t120 * t150 - t125 * t147;
t1 = [qJDD(1) * MDP(1) + t169 * MDP(2) + t170 * MDP(3) + (-0.2e1 * t214 + t222) * MDP(4) + t159 * MDP(5) + (t220 * pkin(1) + (t161 + t183) * qJ(2)) * MDP(6) + (qJDD(1) * t146 - 0.2e1 * t148 * t177) * MDP(7) + 0.2e1 * (-t148 * t181 + t197 * t185) * MDP(8) + (-t148 * t154 + t142) * MDP(9) + (-qJDD(3) * t148 - t151 * t154) * MDP(10) + (t157 * t148 + t163 * t151) * MDP(12) + (-t163 * t148 + t157 * t151) * MDP(13) + (t114 * t150 * t151 + t162 * t131) * MDP(14) + ((t129 * t150 + t131 * t147) * t191 + (-t212 - t115 * t150 + (t129 * t147 - t131 * t150) * qJD(4)) * t151) * MDP(15) + ((t114 - t178) * t148 + (t165 + t193) * t151) * MDP(16) + ((-t115 + t179) * t148 + (-t166 - t194) * t151) * MDP(17) + (t127 * t148 + t138 * t190) * MDP(18) + (-t115 * t199 - g(1) * t124 - g(2) * t122 + t116 * t148 + (t112 * t151 + t129 * t203) * qJD(3) + (t209 + t128 * t138 + (t126 * t151 + (-t138 * t153 - t125) * t148) * qJD(4)) * t150 + t219 * t147) * MDP(19) + (-t114 * t199 + g(1) * t123 - g(2) * t121 + (-t113 * t151 + t131 * t203) * qJD(3) + (-(-qJD(4) * t203 + t128) * t138 - t209 + t173 * t148 - t126 * t187) * t147 + t219 * t150) * MDP(20); qJDD(1) * MDP(4) - t155 * MDP(5) + (-t213 - t220) * MDP(6) + t142 * MDP(12) + (-t150 * MDP(19) + t147 * MDP(20)) * t138 * qJD(1) + (t196 * MDP(13) + (-t115 - t179) * MDP(19) + (-t114 - t178) * MDP(20)) * t151 + (t196 * MDP(12) - qJDD(3) * MDP(13) + (-t166 + t194) * MDP(19) + (-t165 + t193) * MDP(20)) * t148; MDP(9) * t181 - MDP(10) * t182 + qJDD(3) * MDP(11) + (-t160 * t151 + t216) * MDP(12) + (t160 * t148 + t215) * MDP(13) + (t131 * t201 + t212) * MDP(14) + ((t114 - t211) * t150 + (-t115 - t210) * t147) * MDP(15) + ((-t131 * t151 + t148 * t201) * qJD(1) + t166) * MDP(16) + ((t129 * t151 - t138 * t206) * qJD(1) + t165) * MDP(17) - t138 * MDP(18) * t195 + (-pkin(3) * t115 - t112 * t195 - t129 * t208 + t158 * t147 - t218 * t150) * MDP(19) + (-pkin(3) * t114 + t113 * t195 - t131 * t208 + t218 * t147 + t158 * t150) * MDP(20) + (t151 * t148 * MDP(7) - t197 * MDP(8)) * t155; t131 * t129 * MDP(14) + (-t129 ^ 2 + t131 ^ 2) * MDP(15) + (t180 + t211) * MDP(16) + (-t198 + t210) * MDP(17) + t127 * MDP(18) + (-g(1) * t121 - g(2) * t123 + t113 * t138 - t126 * t131 + t116) * MDP(19) + (g(1) * t122 - g(2) * t124 + t112 * t138 + t126 * t129) * MDP(20) + (-MDP(19) * t189 + t168 * MDP(20) + (-MDP(16) * t191 - MDP(17) * t187) * qJD(1)) * t150 + (-qJD(1) * MDP(16) * t187 - t221 * MDP(17) + t168 * MDP(19) + t173 * MDP(20)) * t147;];
tau = t1;
