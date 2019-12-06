% Calculate vector of inverse dynamics joint torques for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PRPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:45
% EndTime: 2019-12-05 15:33:49
% DurationCPUTime: 1.17s
% Computational Cost: add. (645->171), mult. (1293->237), div. (0->0), fcn. (877->10), ass. (0->85)
t220 = qJ(5) + pkin(6);
t164 = sin(qJ(2));
t166 = cos(qJ(2));
t190 = qJD(1) * qJD(2);
t219 = qJDD(1) * t164 + t166 * t190;
t155 = qJ(2) + pkin(8);
t149 = sin(t155);
t150 = cos(t155);
t159 = sin(pkin(7));
t161 = cos(pkin(7));
t179 = g(1) * t161 + g(2) * t159;
t172 = g(3) * t149 + t179 * t150;
t215 = -g(1) * t159 + g(2) * t161;
t218 = qJDD(3) + t215;
t163 = sin(qJ(4));
t156 = t163 ^ 2;
t165 = cos(qJ(4));
t157 = t165 ^ 2;
t216 = (t156 - t157) * MDP(7);
t194 = qJD(1) * t166;
t140 = qJD(2) * pkin(2) + t194;
t160 = cos(pkin(8));
t195 = qJD(1) * t164;
t144 = t160 * t195;
t158 = sin(pkin(8));
t124 = t158 * t140 + t144;
t180 = t220 * qJD(2) + t124;
t114 = t165 * qJD(3) - t180 * t163;
t206 = qJD(4) * pkin(4);
t111 = t114 + t206;
t115 = qJD(3) * t163 + t180 * t165;
t214 = -t111 * t163 + t115 * t165;
t213 = MDP(13) * (t156 + t157);
t208 = g(3) * t150;
t207 = g(3) * t166;
t203 = t149 * t159;
t202 = t149 * t161;
t201 = t150 * t220;
t146 = pkin(2) * t158 + pkin(6);
t200 = qJ(5) + t146;
t199 = qJDD(1) - g(3);
t198 = -t114 + t111;
t193 = qJD(2) * t165;
t192 = qJD(4) * t163;
t191 = qJD(4) * t165;
t188 = qJDD(2) * t163;
t187 = qJDD(2) * t165;
t186 = qJDD(4) * t146;
t152 = t166 * qJDD(1);
t130 = qJDD(2) * pkin(2) - t164 * t190 + t152;
t117 = t158 * t130 + t219 * t160;
t148 = pkin(4) * t165 + pkin(3);
t143 = t158 * t195;
t123 = t140 * t160 - t143;
t181 = qJD(4) * t200;
t116 = t160 * t130 - t219 * t158;
t119 = -qJD(2) * pkin(3) - t123;
t129 = t160 * t194 - t143;
t147 = -pkin(2) * t160 - pkin(3);
t177 = qJD(2) * t147 + t119 + t129;
t175 = t158 * t166 + t160 * t164;
t133 = t158 * t164 - t160 * t166;
t174 = t180 * qJD(4);
t167 = qJD(4) ^ 2;
t137 = qJDD(4) * t165 - t163 * t167;
t136 = qJDD(4) * t163 + t165 * t167;
t173 = t215 * t165;
t113 = qJDD(2) * pkin(6) + t117;
t171 = qJ(5) * qJDD(2) + qJD(2) * qJD(5) + qJD(3) * qJD(4) + t113;
t110 = qJD(2) * pkin(4) * t192 - t148 * qJDD(2) + qJDD(5) - t116;
t170 = -t119 * qJD(2) - t113 + t172;
t127 = t158 * t194 + t144;
t169 = g(1) * t202 + g(2) * t203 + qJD(2) * t127 - t146 * t167 + t116 - t208 + (pkin(3) - t147) * qJDD(2);
t168 = qJD(2) ^ 2;
t151 = t165 * qJDD(3);
t132 = t200 * t165;
t131 = t200 * t163;
t128 = t133 * qJD(2);
t126 = t175 * qJD(2);
t122 = -qJD(5) * t163 - t165 * t181;
t121 = qJD(5) * t165 - t163 * t181;
t118 = -t148 * qJD(2) + qJD(5) - t123;
t109 = (qJDD(3) - t174) * t163 + t171 * t165;
t108 = qJDD(4) * pkin(4) - t171 * t163 - t165 * t174 + t151;
t1 = [t199 * MDP(1) + (qJDD(2) * t166 - t164 * t168) * MDP(3) + (-qJDD(2) * t164 - t166 * t168) * MDP(4) + (-t116 * t133 - t123 * t126 - t124 * t128 - g(3)) * MDP(5) + (t128 * t192 - t133 * t187) * MDP(11) + (t128 * t191 + t133 * t188) * MDP(12) + (t110 * t133 + t118 * t126 - t214 * t128 - g(3)) * MDP(14) + ((-t126 * t165 + t133 * t192) * MDP(11) + (t126 * t163 + t133 * t191) * MDP(12) - t128 * t213) * qJD(2) + (t117 * MDP(5) - t136 * MDP(11) - t137 * MDP(12) + (-t108 * t163 + t109 * t165 - t111 * t191 - t115 * t192) * MDP(14) + qJDD(2) * t213) * t175; t152 * MDP(3) + (t123 * t127 - t124 * t129) * MDP(5) + t136 * MDP(8) + t137 * MDP(9) - t172 * MDP(13) + (t109 * t132 + t115 * t121 - t108 * t131 + t111 * t122 - t110 * pkin(3) - t118 * t127 - g(1) * (-t148 * t202 + t161 * t201) - g(2) * (-t148 * t203 + t159 * t201) - g(3) * (t148 * t150 + t149 * t220)) * MDP(14) + (-g(3) * MDP(3) + t179 * MDP(4)) * t166 + (t179 * MDP(3) - t199 * MDP(4)) * t164 + (MDP(6) * t156 + MDP(2)) * qJDD(2) + (-0.2e1 * qJD(4) * t216 - t129 * t213) * qJD(2) + ((t116 * t160 + t117 * t158 - t207) * MDP(5) + (-t110 * t160 - t207) * MDP(14) + (MDP(14) + MDP(5)) * t164 * t179) * pkin(2) + (t169 * MDP(11) - MDP(12) * t186 + (qJD(2) * t121 + qJDD(2) * t132 + t109) * MDP(13) + (-pkin(4) * t110 - t115 * t129) * MDP(14) + (t177 * MDP(12) + (qJD(2) * t131 - t111) * MDP(13)) * qJD(4)) * t165 + (0.2e1 * MDP(7) * t187 - MDP(11) * t186 - t169 * MDP(12) + (-qJD(2) * t122 + qJDD(2) * t131 - t108) * MDP(13) + t111 * t129 * MDP(14) + (0.2e1 * MDP(6) * t193 + t177 * MDP(11) + (-qJD(2) * t132 - t115) * MDP(13) + pkin(4) * t118 * MDP(14)) * qJD(4)) * t163; t218 * MDP(5) + t137 * MDP(11) - t136 * MDP(12) + (qJD(4) * t214 + t108 * t165 + t109 * t163 + t215) * MDP(14); MDP(8) * t188 + MDP(9) * t187 + qJDD(4) * MDP(10) + (t170 * t163 + t151 + t173) * MDP(11) + (-t218 * t163 + t170 * t165) * MDP(12) + (-pkin(4) * t188 + (t198 - t206) * t193) * MDP(13) + (t198 * t115 + (t108 + t173 + (-t118 * qJD(2) + t172) * t163) * pkin(4)) * MDP(14) + (-t163 * t165 * MDP(6) + t216) * t168; (-qJD(2) * t214 - t179 * t149 + t110 + t208) * MDP(14) - t168 * t213;];
tau = t1;
