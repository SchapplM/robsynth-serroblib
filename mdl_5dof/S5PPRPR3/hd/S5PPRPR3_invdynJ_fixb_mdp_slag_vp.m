% Calculate vector of inverse dynamics joint torques for
% S5PPRPR3
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:26
% EndTime: 2019-12-05 15:05:28
% DurationCPUTime: 0.63s
% Computational Cost: add. (396->129), mult. (866->196), div. (0->0), fcn. (715->12), ass. (0->76)
t149 = sin(pkin(8));
t150 = sin(pkin(7));
t153 = cos(pkin(7));
t200 = (g(1) * t153 + g(2) * t150) * t149;
t155 = sin(qJ(3));
t157 = cos(qJ(3));
t181 = qJD(1) * t149;
t199 = t157 * qJD(2) - t155 * t181;
t148 = sin(pkin(9));
t151 = cos(pkin(9));
t129 = t148 * t155 - t151 * t157;
t195 = g(3) * t149;
t132 = qJD(2) * t155 + t157 * t181;
t194 = t132 * t148;
t152 = cos(pkin(8));
t192 = t150 * t152;
t191 = t150 * t155;
t190 = t150 * t157;
t189 = t151 * t132;
t187 = t152 * t153;
t186 = t153 * t155;
t185 = t153 * t157;
t154 = sin(qJ(5));
t156 = cos(qJ(5));
t184 = t154 * t156;
t183 = qJDD(1) - g(3);
t142 = t157 * qJDD(2);
t177 = qJDD(1) * t149;
t114 = qJDD(3) * pkin(3) - t132 * qJD(3) - t155 * t177 + t142;
t175 = t155 * qJDD(2);
t115 = qJD(3) * t199 + t157 * t177 + t175;
t105 = t148 * t114 + t151 * t115;
t146 = t154 ^ 2;
t182 = -t156 ^ 2 + t146;
t180 = qJD(3) * t149;
t124 = t129 * t149;
t179 = qJD(5) * t124;
t178 = qJD(3) * qJD(5);
t176 = t124 * qJDD(5);
t173 = -g(1) * t150 + g(2) * t153;
t172 = t149 * t183;
t104 = t114 * t151 - t115 * t148;
t126 = qJD(3) * pkin(3) + t199;
t108 = t126 * t151 - t194;
t130 = t148 * t157 + t151 * t155;
t159 = qJD(3) ^ 2;
t171 = qJDD(3) * t157 - t155 * t159;
t170 = -qJDD(3) * t155 - t157 * t159;
t117 = t130 * t180;
t123 = t130 * t149;
t168 = qJD(3) * t123 + qJD(5) * t152 + t117;
t127 = t130 * qJD(3);
t158 = qJD(5) ^ 2;
t167 = qJD(3) * t127 + qJDD(3) * t129 + t130 * t158;
t116 = t129 * t180;
t166 = qJD(3) * t116 - qJDD(3) * t123 - qJDD(5) * t152;
t128 = qJD(3) * t129;
t165 = qJD(5) * t128 - qJDD(5) * t130 + t129 * t178;
t106 = -qJD(3) * pkin(4) - t108;
t112 = t151 * t199 - t194;
t138 = pkin(3) * t148 + pkin(6);
t139 = -pkin(3) * t151 - pkin(4);
t164 = -qJDD(5) * t138 + (qJD(3) * t139 + t106 + t112) * qJD(5);
t136 = -qJDD(1) * t152 + qJDD(4);
t163 = -g(3) * t152 - t136 + t200;
t145 = qJ(3) + pkin(9);
t140 = sin(t145);
t141 = cos(t145);
t162 = -qJDD(3) * pkin(6) + g(1) * (t140 * t150 + t141 * t187) + g(2) * (-t140 * t153 + t141 * t192) - qJD(3) * t106 + t141 * t195 - t105;
t161 = -g(1) * (-t152 * t186 + t190) - g(2) * (-t152 * t191 - t185);
t111 = t148 * t199 + t189;
t160 = -g(1) * (-t140 * t187 + t141 * t150) - g(2) * (-t140 * t192 - t141 * t153) + qJD(3) * t111 - t138 * t158 + t140 * t195 + t104 + (pkin(4) - t139) * qJDD(3);
t135 = qJDD(5) * t156 - t154 * t158;
t134 = qJDD(5) * t154 + t156 * t158;
t109 = t148 * t126 + t189;
t1 = [t183 * MDP(1) + (-g(3) + (t149 ^ 2 + t152 ^ 2) * qJDD(1)) * MDP(2) + t170 * t149 * MDP(4) - t171 * t149 * MDP(5) + (-t104 * t123 - t105 * t124 + t108 * t116 - t109 * t117 - t136 * t152 - g(3)) * MDP(6) + (t154 * t176 + t166 * t156 + (t168 * t154 + t156 * t179) * qJD(5)) * MDP(12) + (t156 * t176 - t166 * t154 + (-t154 * t179 + t168 * t156) * qJD(5)) * MDP(13); (qJDD(2) + t173) * MDP(2) + t171 * MDP(4) + t170 * MDP(5) + (-t104 * t129 + t105 * t130 - t108 * t127 - t109 * t128 + t173) * MDP(6) + (-t167 * MDP(12) + t165 * MDP(13)) * t156 + (t165 * MDP(12) + t167 * MDP(13)) * t154; qJDD(3) * MDP(3) + (-t155 * t172 + t142 + t161) * MDP(4) + (-t175 - g(1) * (-t152 * t185 - t191) - g(2) * (-t152 * t190 + t186) - t157 * t172) * MDP(5) + (t108 * t111 - t109 * t112 + (t104 * t151 + t105 * t148 + t155 * t195 + t161) * pkin(3)) * MDP(6) + (qJDD(3) * t146 + 0.2e1 * t178 * t184) * MDP(7) + 0.2e1 * (qJDD(3) * t184 - t182 * t178) * MDP(8) + t134 * MDP(9) + t135 * MDP(10) + (t164 * t154 + t160 * t156) * MDP(12) + (-t160 * t154 + t164 * t156) * MDP(13); (-t183 * t152 + qJDD(4) - t200) * MDP(6) + t135 * MDP(12) - t134 * MDP(13); qJDD(5) * MDP(11) + t182 * MDP(8) * t159 + (qJDD(3) * MDP(10) - t163 * MDP(12) + t162 * MDP(13)) * t156 + (-t159 * t156 * MDP(7) + t162 * MDP(12) + t163 * MDP(13) + qJDD(3) * MDP(9)) * t154;];
tau = t1;
