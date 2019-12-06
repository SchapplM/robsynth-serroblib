% Calculate vector of inverse dynamics joint torques for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:18
% EndTime: 2019-12-05 15:07:21
% DurationCPUTime: 1.00s
% Computational Cost: add. (545->143), mult. (1169->197), div. (0->0), fcn. (830->10), ass. (0->71)
t139 = sin(pkin(8));
t141 = cos(pkin(8));
t145 = sin(qJ(3));
t147 = cos(qJ(3));
t207 = -t139 * t145 + t141 * t147;
t115 = t207 * qJD(1);
t160 = t139 * t147 + t141 * t145;
t118 = t160 * qJD(3);
t192 = qJ(5) + pkin(6);
t136 = pkin(8) + qJ(3);
t131 = sin(t136);
t132 = cos(t136);
t140 = sin(pkin(7));
t142 = cos(pkin(7));
t164 = g(1) * t142 + g(2) * t140;
t206 = -g(3) * t132 + t164 * t131;
t153 = g(3) * t131 + t164 * t132;
t200 = -g(1) * t140 + g(2) * t142;
t205 = qJDD(2) + t200;
t144 = sin(qJ(4));
t137 = t144 ^ 2;
t146 = cos(qJ(4));
t138 = t146 ^ 2;
t202 = (t137 - t138) * MDP(7);
t167 = (t137 + t138) * MDP(13);
t201 = MDP(5) - t167;
t172 = qJDD(1) * t147;
t173 = qJDD(1) * t145;
t156 = qJD(1) * t118 + t139 * t173 - t141 * t172;
t199 = -t156 + t206;
t116 = t160 * qJD(1);
t166 = t192 * qJD(3) + t116;
t108 = t146 * qJD(2) - t166 * t144;
t190 = qJD(4) * pkin(4);
t107 = t108 + t190;
t109 = qJD(2) * t144 + t166 * t146;
t198 = -t107 * t144 + t109 * t146;
t191 = qJD(3) * pkin(3);
t189 = pkin(6) * qJDD(4);
t181 = -t108 + t107;
t178 = MDP(14) * t115;
t177 = qJD(3) * t146;
t176 = qJD(4) * t144;
t175 = qJD(4) * t146;
t171 = qJDD(3) * t144;
t170 = qJDD(3) * t146;
t130 = pkin(4) * t146 + pkin(3);
t168 = qJD(4) * t192;
t111 = -t115 - t191;
t165 = t111 + t115 - t191;
t158 = t166 * qJD(4);
t148 = qJD(4) ^ 2;
t122 = qJDD(4) * t146 - t144 * t148;
t121 = qJDD(4) * t144 + t146 * t148;
t157 = t200 * t146;
t117 = t207 * qJD(3);
t105 = qJDD(3) * pkin(6) + qJD(1) * t117 + t160 * qJDD(1);
t152 = qJ(5) * qJDD(3) + qJD(2) * qJD(4) + qJD(3) * qJD(5) + t105;
t151 = -t111 * qJD(3) - t105 + t153;
t104 = qJD(3) * pkin(4) * t176 - t130 * qJDD(3) + qJDD(5) + t156;
t150 = 0.2e1 * qJDD(3) * pkin(3) - pkin(6) * t148 + qJD(3) * t116 + t199;
t149 = qJD(3) ^ 2;
t133 = t146 * qJDD(2);
t124 = t192 * t146;
t123 = t192 * t144;
t114 = -qJD(5) * t144 - t146 * t168;
t113 = qJD(5) * t146 - t144 * t168;
t110 = -t130 * qJD(3) + qJD(5) - t115;
t103 = (qJDD(2) - t158) * t144 + t152 * t146;
t102 = qJDD(4) * pkin(4) - t152 * t144 - t146 * t158 + t133;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t139 ^ 2 + t141 ^ 2) * qJDD(1)) * MDP(2) + t207 * qJDD(3) * MDP(4) + (-t117 * t176 + t170 * t207) * MDP(11) + (-t117 * t175 - t171 * t207) * MDP(12) + (-t104 * t207 + t110 * t118 + t198 * t117 - g(3)) * MDP(14) + (-t118 * MDP(4) + (-t118 * t146 - t176 * t207) * MDP(11) + (t118 * t144 - t175 * t207) * MDP(12) - t201 * t117) * qJD(3) + (-t121 * MDP(11) - t122 * MDP(12) + (-t102 * t144 + t103 * t146 - t107 * t175 - t109 * t176) * MDP(14) - t201 * qJDD(3)) * t160; t205 * MDP(2) + t122 * MDP(11) - t121 * MDP(12) + (qJD(4) * t198 + t102 * t146 + t103 * t144 + t200) * MDP(14); t199 * MDP(4) + (-t139 * t172 - t141 * t173 + t153) * MDP(5) + t121 * MDP(8) + t122 * MDP(9) - t153 * MDP(13) + (t103 * t124 + t109 * t113 - t102 * t123 + t107 * t114 - t104 * t130 - t110 * t116 - g(3) * (t130 * t132 + t131 * t192) + t164 * (t130 * t131 - t132 * t192)) * MDP(14) + (t137 * MDP(6) + MDP(3)) * qJDD(3) + (t116 * MDP(4) - 0.2e1 * qJD(4) * t202 + (-MDP(5) + t201) * t115) * qJD(3) + (t150 * MDP(11) - MDP(12) * t189 + (qJD(3) * t113 + qJDD(3) * t124 + t103) * MDP(13) - t109 * t178 + (t165 * MDP(12) + (qJD(3) * t123 - t107) * MDP(13)) * qJD(4)) * t146 + (0.2e1 * MDP(7) * t170 - MDP(11) * t189 - t150 * MDP(12) + (-qJD(3) * t114 + qJDD(3) * t123 - t102) * MDP(13) + t107 * t178 + (0.2e1 * MDP(6) * t177 + t165 * MDP(11) + (-qJD(3) * t124 - t109) * MDP(13) + pkin(4) * t110 * MDP(14)) * qJD(4)) * t144; MDP(8) * t171 + MDP(9) * t170 + qJDD(4) * MDP(10) + (t151 * t144 + t133 + t157) * MDP(11) + (-t205 * t144 + t151 * t146) * MDP(12) + (-pkin(4) * t171 + (t181 - t190) * t177) * MDP(13) + (t181 * t109 + (t102 + t157 + (-t110 * qJD(3) + t153) * t144) * pkin(4)) * MDP(14) + (-t144 * t146 * MDP(6) + t202) * t149; (-qJD(3) * t198 + t104 - t206) * MDP(14) - t149 * t167;];
tau = t1;
