% Calculate vector of inverse dynamics joint torques for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:25
% EndTime: 2019-12-05 15:03:27
% DurationCPUTime: 0.59s
% Computational Cost: add. (315->107), mult. (613->134), div. (0->0), fcn. (448->10), ass. (0->63)
t127 = sin(pkin(8));
t129 = cos(pkin(8));
t132 = sin(qJ(3));
t134 = cos(qJ(3));
t107 = t127 * t134 + t129 * t132;
t103 = t107 * qJD(1);
t167 = qJD(3) * qJ(4);
t101 = t103 + t167;
t163 = qJD(1) * qJD(3);
t154 = t132 * t163;
t160 = qJDD(1) * t134;
t181 = qJDD(1) * t132 + t134 * t163;
t145 = (-t154 + t160) * t129 - t181 * t127;
t143 = qJDD(4) - t145;
t122 = pkin(8) + qJ(3);
t120 = sin(t122);
t121 = cos(t122);
t128 = sin(pkin(7));
t130 = cos(pkin(7));
t152 = g(1) * t130 + g(2) * t128;
t179 = g(3) * t121 - t152 * t120;
t139 = t143 + t179;
t183 = -t101 * qJD(3) + t139;
t172 = t127 * t132;
t155 = qJD(1) * t172;
t171 = t129 * t134;
t102 = -qJD(1) * t171 + t155;
t142 = -g(3) * t120 - t121 * t152;
t156 = t127 * t160 + t129 * t181;
t182 = -t142 + (-t102 + t155) * qJD(3) - t156;
t106 = -t171 + t172;
t135 = -pkin(3) - pkin(6);
t177 = (t101 - t103 + t167) * qJD(5) + qJDD(5) * t135;
t173 = qJDD(3) * pkin(3);
t131 = sin(qJ(5));
t133 = cos(qJ(5));
t170 = t131 * t133;
t126 = t133 ^ 2;
t169 = t131 ^ 2 - t126;
t136 = qJD(5) ^ 2;
t137 = qJD(3) ^ 2;
t168 = -t136 - t137;
t166 = qJD(3) * t103;
t164 = qJD(4) + t102;
t162 = qJD(3) * qJD(5);
t123 = qJDD(3) * qJ(4);
t159 = qJDD(5) * t131;
t158 = qJDD(5) * t133;
t104 = t106 * qJD(3);
t148 = qJD(3) * t104 - qJDD(3) * t107;
t146 = g(1) * t128 - g(2) * t130 - qJDD(2);
t144 = -t106 * t136 - t148;
t105 = qJD(3) * t107;
t141 = qJD(5) * t105 + qJDD(5) * t106 + t107 * t162;
t140 = -t135 * qJDD(3) - t183;
t124 = qJD(3) * qJD(4);
t97 = t127 * t154 - t123 - t124 - t156;
t138 = qJD(3) * t164 - t135 * t136 + t123 + t142 - t97;
t110 = -t131 * t136 + t158;
t109 = -t133 * t136 - t159;
t100 = -qJD(3) * pkin(3) + t164;
t98 = t143 - t173;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t127 ^ 2 + t129 ^ 2) * qJDD(1)) * MDP(2) + (t100 * t105 - t101 * t104 + t106 * t98 - t107 * t97 - g(3)) * MDP(8) + (MDP(14) * t141 + MDP(15) * t144) * t133 + (MDP(14) * t144 - MDP(15) * t141) * t131 + (MDP(5) - MDP(7)) * t148 + (-MDP(4) + MDP(6)) * (qJD(3) * t105 + qJDD(3) * t106); MDP(14) * t109 - MDP(15) * t110 - (MDP(2) + MDP(8)) * t146; qJDD(3) * MDP(3) + (t145 + t166 - t179) * MDP(4) + t182 * MDP(5) + (t139 - t166 - 0.2e1 * t173) * MDP(6) + (0.2e1 * t123 + 0.2e1 * t124 - t182) * MDP(7) + (-t97 * qJ(4) - t98 * pkin(3) - t100 * t103 - g(3) * (pkin(3) * t121 + qJ(4) * t120) + t164 * t101 + t152 * (pkin(3) * t120 - qJ(4) * t121)) * MDP(8) + (qJDD(3) * t126 - 0.2e1 * t162 * t170) * MDP(9) + 0.2e1 * (-qJDD(3) * t170 + t162 * t169) * MDP(10) + t110 * MDP(11) + t109 * MDP(12) + (t138 * t131 + t133 * t177) * MDP(14) + (-t131 * t177 + t138 * t133) * MDP(15); qJDD(3) * MDP(6) - t137 * MDP(7) + (-t173 + t183) * MDP(8) + (t131 * t168 + t158) * MDP(14) + (t133 * t168 - t159) * MDP(15); qJDD(5) * MDP(13) - t169 * MDP(10) * t137 + (qJDD(3) * MDP(11) - MDP(14) * t140 + MDP(15) * t146) * t133 + (t133 * t137 * MDP(9) - qJDD(3) * MDP(12) + MDP(14) * t146 + MDP(15) * t140) * t131;];
tau = t1;
