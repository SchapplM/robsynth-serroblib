% Calculate vector of inverse dynamics joint torques for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:46
% EndTime: 2019-12-31 16:40:48
% DurationCPUTime: 0.63s
% Computational Cost: add. (324->159), mult. (735->212), div. (0->0), fcn. (487->6), ass. (0->82)
t147 = cos(qJ(1));
t139 = g(2) * t147;
t145 = sin(qJ(1));
t190 = g(1) * t145;
t162 = -t139 + t190;
t143 = cos(pkin(6));
t191 = pkin(2) * t143;
t189 = -pkin(5) + qJ(2);
t188 = qJDD(1) * pkin(1);
t142 = sin(pkin(6));
t140 = t142 ^ 2;
t149 = qJD(1) ^ 2;
t187 = t140 * t149;
t146 = cos(qJ(4));
t186 = t142 * t146;
t144 = sin(qJ(4));
t185 = t143 * t144;
t184 = t143 * t145;
t141 = t143 ^ 2;
t176 = (qJD(1) * qJD(2));
t169 = 2 * t176;
t129 = t141 * t169;
t170 = qJDD(1) * qJ(2) ^ 2;
t183 = qJ(2) * t129 + t141 * t170;
t182 = t147 * pkin(1) + t145 * qJ(2);
t181 = qJD(1) * t142;
t180 = qJD(1) * t144;
t179 = qJD(1) * t146;
t178 = qJD(3) * t142;
t177 = qJ(2) * qJDD(1);
t175 = qJD(1) * qJD(3);
t164 = qJ(3) * t142 + pkin(1);
t118 = -t164 - t191;
t174 = qJDD(1) * t118;
t173 = qJDD(1) * t142;
t172 = qJDD(1) * t144;
t171 = qJDD(1) * t146;
t168 = t142 * t179;
t167 = qJD(4) * t186;
t166 = t143 * t180;
t165 = t142 * t175;
t113 = qJ(2) * t173 + t142 * t176 + qJDD(3);
t163 = -pkin(1) * t145 + t147 * qJ(2);
t161 = g(1) * t147 + g(2) * t145;
t125 = t142 * t171;
t160 = -t143 * t172 + t125;
t159 = -qJDD(2) + t165;
t120 = t189 * t142;
t121 = t189 * t143;
t158 = t120 * t146 - t121 * t144;
t157 = t120 * t144 + t121 * t146;
t116 = -t185 + t186;
t115 = t142 * t144 + t143 * t146;
t135 = qJDD(2) - t188;
t156 = -t135 + t188 - t139;
t98 = -t159 + t174;
t155 = -t174 - t98 - t139;
t154 = 0.2e1 * t141 * t177 + t129 - t161;
t107 = t115 * qJD(4);
t112 = (pkin(2) + pkin(3)) * t143 + t164;
t153 = (t176 + t177) * t140;
t152 = -t142 * t180 - t143 * t179;
t151 = -qJ(2) * t141 * t149 + qJDD(2) - t162;
t148 = qJD(4) ^ 2;
t130 = g(1) * t184;
t123 = qJ(2) * t181 + qJD(3);
t122 = qJD(4) * t166;
t111 = -t166 + t168;
t109 = t115 * qJD(1);
t108 = -qJD(4) * t185 + t167;
t106 = t118 * qJD(1) + qJD(2);
t105 = t115 * t147;
t104 = t116 * t147;
t103 = t115 * t145;
t102 = t116 * t145;
t101 = (t189 * qJDD(1) + t176) * t143;
t100 = -pkin(5) * t173 + t113;
t99 = t112 * qJD(1) - qJD(2);
t97 = t112 * qJDD(1) + t159;
t96 = qJD(1) * t167 + t115 * qJDD(1) - t122;
t95 = -qJD(1) * t107 + t160;
t1 = [qJDD(1) * MDP(1) + t162 * MDP(2) + t161 * MDP(3) + (t156 * t143 + t130) * MDP(4) + (-t156 - t190) * t142 * MDP(5) + (0.2e1 * t153 + t154) * MDP(6) + (-t135 * pkin(1) - g(1) * t163 - g(2) * t182 + (qJ(2) * t169 + t170) * t140 + t183) * MDP(7) + (t130 + (t155 + t165) * t143) * MDP(8) + (t113 * t142 + t153 + t154) * MDP(9) + (t140 * t175 + (t155 + t190) * t142) * MDP(10) + (t98 * t118 - g(1) * (-pkin(2) * t184 + t163) - g(2) * (t147 * t191 + t182) + (t113 * qJ(2) + t162 * qJ(3) + t123 * qJD(2) - t106 * qJD(3)) * t142 + t183) * MDP(11) + (-t107 * t111 + t116 * t95) * MDP(12) + (t107 * t109 - t108 * t111 - t115 * t95 - t116 * t96) * MDP(13) + (-qJD(4) * t107 + qJDD(4) * t116) * MDP(14) + (-qJD(4) * t108 - qJDD(4) * t115) * MDP(15) + (t109 * t178 + t112 * t96 + t97 * t115 + t99 * t108 + t158 * qJDD(4) + g(1) * t103 - g(2) * t105 + (t116 * qJD(2) - t157 * qJD(4)) * qJD(4)) * MDP(17) + (t111 * t178 + t112 * t95 + t97 * t116 - t99 * t107 - t157 * qJDD(4) + g(1) * t102 - g(2) * t104 + (-t115 * qJD(2) - t158 * qJD(4)) * qJD(4)) * MDP(18); (-qJ(2) * t187 + t151) * MDP(7) + (-t123 * t181 + t151 - t165) * MDP(11) + t122 * MDP(17) - t125 * MDP(18) + (MDP(6) + MDP(9)) * (-t140 - t141) * t149 + ((-t111 - t168) * MDP(17) + (t109 - t152) * MDP(18)) * qJD(4) + ((-MDP(7) - MDP(11)) * pkin(1) + (-MDP(11) * qJ(3) - MDP(17) * t144 - MDP(10) + MDP(5)) * t142 + (-MDP(11) * pkin(2) - MDP(17) * t146 + MDP(18) * t144 - MDP(4) - MDP(8)) * t143) * qJDD(1); -MDP(10) * t187 + (g(3) * t143 + t113) * MDP(11) + (t146 * qJDD(4) - t144 * t148) * MDP(17) + (-t144 * qJDD(4) - t146 * t148) * MDP(18) + (-t149 * t143 * MDP(8) + qJDD(1) * MDP(9) - t161 * MDP(11) + (t106 * MDP(11) - t109 * MDP(17) - t111 * MDP(18)) * qJD(1)) * t142; t111 * t109 * MDP(12) + (-t109 ^ 2 + t111 ^ 2) * MDP(13) + t160 * MDP(14) + (-t142 * t172 - t143 * t171 + t122) * MDP(15) + qJDD(4) * MDP(16) + (-g(1) * t104 - g(2) * t102 + g(3) * t115 + t146 * t100 - t144 * t101 - t99 * t111) * MDP(17) + (g(1) * t105 + g(2) * t103 + g(3) * t116 - t144 * t100 - t146 * t101 + t99 * t109) * MDP(18) + ((t109 + t152) * MDP(14) + (t111 - t168) * MDP(15)) * qJD(4);];
tau = t1;
