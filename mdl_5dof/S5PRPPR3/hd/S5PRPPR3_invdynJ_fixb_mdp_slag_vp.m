% Calculate vector of inverse dynamics joint torques for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:59
% EndTime: 2019-12-05 15:27:01
% DurationCPUTime: 0.68s
% Computational Cost: add. (376->130), mult. (701->164), div. (0->0), fcn. (493->10), ass. (0->68)
t124 = qJ(2) + pkin(8);
t121 = sin(t124);
t122 = cos(t124);
t134 = cos(qJ(2));
t123 = t134 * qJDD(1);
t132 = sin(qJ(2));
t158 = qJD(1) * qJD(2);
t106 = qJDD(2) * pkin(2) - t132 * t158 + t123;
t127 = sin(pkin(8));
t129 = cos(pkin(8));
t153 = t132 * qJDD(1);
t177 = t134 * t158 + t153;
t95 = t129 * t106 - t177 * t127;
t146 = qJDD(4) - t95;
t128 = sin(pkin(7));
t130 = cos(pkin(7));
t150 = g(1) * t130 + g(2) * t128;
t138 = g(3) * t122 - t150 * t121 + t146;
t160 = t134 * qJD(1);
t113 = qJD(2) * pkin(2) + t160;
t161 = qJD(1) * t132;
t101 = t127 * t113 + t129 * t161;
t99 = qJD(2) * qJ(4) + t101;
t178 = -t99 * qJD(2) + t138;
t174 = -t127 * t132 + t129 * t134;
t108 = t127 * t134 + t129 * t132;
t103 = t108 * qJD(1);
t120 = -t129 * pkin(2) - pkin(3);
t117 = -pkin(6) + t120;
t118 = t127 * pkin(2) + qJ(4);
t173 = qJDD(5) * t117 + (qJD(2) * t118 - t103 + t99) * qJD(5);
t170 = qJDD(2) * pkin(3);
t169 = t127 * t106;
t131 = sin(qJ(5));
t133 = cos(qJ(5));
t166 = t131 * t133;
t164 = qJDD(1) - g(3);
t126 = t133 ^ 2;
t163 = t131 ^ 2 - t126;
t135 = qJD(5) ^ 2;
t136 = qJD(2) ^ 2;
t162 = -t135 - t136;
t115 = t127 * t161;
t152 = t129 * t160;
t105 = -t115 + t152;
t159 = qJD(4) - t105;
t157 = qJD(2) * qJD(5);
t155 = qJDD(5) * t131;
t154 = qJDD(5) * t133;
t100 = t129 * t113 - t115;
t104 = t174 * qJD(2);
t148 = t104 * qJD(2) + t108 * qJDD(2);
t147 = g(1) * t128 - g(2) * t130 - qJDD(3);
t145 = t129 * t153 + t169;
t143 = t135 * t174 + t148;
t142 = -g(3) * t121 - t150 * t122;
t141 = -g(3) * t134 + t150 * t132;
t102 = qJD(2) * t108;
t140 = qJD(5) * t102 - qJDD(5) * t174 + t108 * t157;
t139 = -(-pkin(3) - pkin(6)) * qJDD(2) - t178;
t93 = qJDD(2) * qJ(4) + (qJD(4) + t152) * qJD(2) + t145;
t137 = t159 * qJD(2) + t118 * qJDD(2) - t117 * t135 + t142 + t93;
t111 = -t135 * t131 + t154;
t110 = -t135 * t133 - t155;
t98 = -qJD(2) * pkin(3) + qJD(4) - t100;
t96 = t177 * t129 + t169;
t94 = t146 - t170;
t1 = [t164 * MDP(1) + (t134 * qJDD(2) - t136 * t132) * MDP(3) + (-qJDD(2) * t132 - t136 * t134) * MDP(4) + (-t100 * t102 + t101 * t104 + t96 * t108 + t174 * t95 - g(3)) * MDP(5) + (t102 * qJD(2) - qJDD(2) * t174) * MDP(6) + t148 * MDP(7) + (t98 * t102 + t99 * t104 + t93 * t108 - t174 * t94 - g(3)) * MDP(8) + (t140 * MDP(14) + t143 * MDP(15)) * t133 + (t143 * MDP(14) - t140 * MDP(15)) * t131; qJDD(2) * MDP(2) + (t123 + t141) * MDP(3) + (-t164 * t132 + t150 * t134) * MDP(4) + (t100 * t103 - t101 * t105 + (t127 * t96 + t129 * t95 + t141) * pkin(2)) * MDP(5) + (-t103 * qJD(2) + (-pkin(3) + t120) * qJDD(2) + t138) * MDP(6) + ((qJ(4) + t118) * qJDD(2) + (0.2e1 * qJD(4) - t105 + t152) * qJD(2) + t142 + t145) * MDP(7) + (t93 * t118 + t94 * t120 - t98 * t103 - g(3) * (t134 * pkin(2) + t122 * pkin(3) + t121 * qJ(4)) + t159 * t99 + t150 * (pkin(2) * t132 + pkin(3) * t121 - qJ(4) * t122)) * MDP(8) + (t126 * qJDD(2) - 0.2e1 * t157 * t166) * MDP(9) + 0.2e1 * (-qJDD(2) * t166 + t163 * t157) * MDP(10) + t111 * MDP(11) + t110 * MDP(12) + (t137 * t131 + t173 * t133) * MDP(14) + (-t173 * t131 + t137 * t133) * MDP(15); t110 * MDP(14) - t111 * MDP(15) - (MDP(5) + MDP(8)) * t147; qJDD(2) * MDP(6) - t136 * MDP(7) + (-t170 + t178) * MDP(8) + (t162 * t131 + t154) * MDP(14) + (t162 * t133 - t155) * MDP(15); qJDD(5) * MDP(13) - t163 * MDP(10) * t136 + (qJDD(2) * MDP(11) - t139 * MDP(14) + t147 * MDP(15)) * t133 + (t133 * t136 * MDP(9) - qJDD(2) * MDP(12) + t147 * MDP(14) + t139 * MDP(15)) * t131;];
tau = t1;
