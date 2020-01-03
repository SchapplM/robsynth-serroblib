% Calculate vector of inverse dynamics joint torques for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:26
% EndTime: 2019-12-31 17:32:28
% DurationCPUTime: 0.94s
% Computational Cost: add. (407->145), mult. (851->190), div. (0->0), fcn. (647->10), ass. (0->82)
t149 = sin(pkin(8));
t150 = cos(pkin(8));
t192 = t149 ^ 2 + t150 ^ 2;
t195 = t150 * MDP(6);
t213 = -MDP(7) * t149 + t195;
t151 = sin(qJ(5));
t153 = cos(qJ(5));
t124 = t149 * t153 + t150 * t151;
t211 = t124 * qJD(5);
t212 = qJD(3) * t211;
t152 = sin(qJ(3));
t191 = qJD(2) * t152;
t169 = t192 * (qJD(3) * qJ(4) + t191);
t154 = cos(qJ(3));
t203 = sin(pkin(7));
t204 = cos(pkin(7));
t121 = -t203 * t152 - t204 * t154;
t123 = t204 * t152 - t203 * t154;
t173 = g(1) * t121 + g(2) * t123;
t190 = qJD(2) * t154;
t175 = qJD(4) - t190;
t141 = qJD(3) * t191;
t210 = qJDD(4) + t141;
t209 = -t169 * t154 - (-qJD(3) * pkin(3) + t175) * t152;
t208 = qJD(5) ^ 2;
t205 = pkin(6) + qJ(4);
t201 = pkin(6) * qJDD(3);
t200 = qJDD(3) * pkin(3);
t196 = t149 * t151;
t194 = t153 * t150;
t155 = qJD(3) ^ 2;
t193 = t154 * t155;
t189 = qJD(3) * t151;
t188 = qJD(3) * t153;
t187 = qJDD(1) * t150;
t186 = qJDD(2) * t152;
t185 = qJDD(2) * t154;
t184 = qJDD(3) * t151;
t183 = qJDD(3) * t153;
t122 = -t194 + t196;
t182 = qJDD(5) * t122;
t181 = qJDD(5) * t124;
t177 = t150 * t188;
t180 = qJD(5) * t177 + t149 * t183 + t150 * t184;
t140 = -pkin(4) * t150 - pkin(3);
t179 = t149 * t189;
t178 = qJD(5) * t196;
t176 = t192 * qJDD(3);
t111 = qJDD(3) * qJ(4) + t186 + (qJD(4) + t190) * qJD(3);
t106 = -qJDD(1) * t149 + t150 * t111;
t174 = g(1) * t123 - g(2) * t121;
t172 = -g(1) * t203 + g(2) * t204;
t133 = t150 * t183;
t171 = -t149 * t184 + t133;
t105 = -t111 * t149 - t187;
t170 = -t105 * t149 + t106 * t150;
t126 = t205 * t149;
t127 = t205 * t150;
t168 = -t126 * t153 - t127 * t151;
t167 = -t126 * t151 + t127 * t153;
t166 = -t185 + t210;
t118 = -qJD(5) * t194 + t178;
t164 = -qJD(5) * t118 + t181;
t163 = qJD(5) * t211 + t182;
t114 = t166 - t200;
t162 = -t114 + t174;
t160 = t122 * qJD(5);
t159 = -t149 * t188 - t150 * t189;
t158 = t174 + t185;
t157 = t170 + t173;
t147 = pkin(8) + qJ(5);
t143 = cos(t147);
t142 = sin(t147);
t120 = t140 * qJD(3) + t175;
t117 = t124 * qJD(3);
t115 = -t177 + t179;
t109 = t140 * qJDD(3) + t166;
t104 = t150 * t201 + t106;
t103 = -t187 + (-t111 - t201) * t149;
t102 = -t171 + t212;
t101 = -qJD(3) * t178 + t180;
t1 = [(-t105 * t150 - t106 * t149 - g(3)) * MDP(9) + t163 * MDP(15) + t164 * MDP(16) + (MDP(1) + MDP(2)) * (qJDD(1) - g(3)); (qJDD(2) + t172) * MDP(2) + (-qJDD(3) * t152 - t193) * MDP(5) + (t152 * t176 + t192 * t193) * MDP(8) + (-t209 * qJD(3) - t114 * t154 + t170 * t152 + t172) * MDP(9) + ((-t102 - t212) * t154 + (qJD(3) * t115 + t122 * t208 - t181) * t152) * MDP(15) + ((qJD(3) * t160 - t101) * t154 + (qJD(3) * t117 + t124 * t208 + t182) * t152) * MDP(16) + (MDP(4) + t213) * (qJDD(3) * t154 - t152 * t155); qJDD(3) * MDP(3) + t158 * MDP(4) + (-t173 - t186) * MDP(5) + (t175 * qJD(3) * t192 + qJ(4) * t176 + t157) * MDP(8) + (t162 * pkin(3) + t157 * qJ(4) + t209 * qJD(2) + t169 * qJD(4)) * MDP(9) + (t101 * t124 - t117 * t118) * MDP(10) + (-t101 * t122 - t102 * t124 + t115 * t118 - t117 * t211) * MDP(11) + t164 * MDP(12) - t163 * MDP(13) + (t140 * t102 + t109 * t122 + t120 * t211 + t168 * qJDD(5) + t174 * t143 + (-t124 * qJD(4) - t167 * qJD(5)) * qJD(5) + (-t152 * t115 + t154 * t211) * qJD(2)) * MDP(15) + (t140 * t101 + t109 * t124 - t120 * t118 - t167 * qJDD(5) - t174 * t142 + (t122 * qJD(4) - t168 * qJD(5)) * qJD(5) + (-t152 * t117 - t154 * t160) * qJD(2)) * MDP(16) + t213 * (t141 + t162 + t200); (-t169 * qJD(3) - t158 + t210) * MDP(9) - t133 * MDP(15) + t180 * MDP(16) - t192 * MDP(8) * t155 + (-t195 - pkin(3) * MDP(9) + (t151 * MDP(15) + MDP(7)) * t149) * qJDD(3) + ((t117 - t159) * MDP(15) + (-t115 - t179) * MDP(16)) * qJD(5); t117 * t115 * MDP(10) + (-t115 ^ 2 + t117 ^ 2) * MDP(11) + t180 * MDP(12) + t171 * MDP(13) + qJDD(5) * MDP(14) + (g(3) * t143 + t153 * t103 - t151 * t104 - t120 * t117 - t173 * t142) * MDP(15) + (-g(3) * t142 - t151 * t103 - t153 * t104 + t120 * t115 - t173 * t143) * MDP(16) + ((t115 - t179) * MDP(12) + (t117 + t159) * MDP(13)) * qJD(5);];
tau = t1;
