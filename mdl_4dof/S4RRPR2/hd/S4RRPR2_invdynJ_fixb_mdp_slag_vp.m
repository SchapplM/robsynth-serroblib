% Calculate vector of inverse dynamics joint torques for
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:36
% EndTime: 2019-07-18 18:16:37
% DurationCPUTime: 0.45s
% Computational Cost: add. (466->129), mult. (513->160), div. (0->0), fcn. (246->8), ass. (0->62)
t123 = qJ(1) + qJ(2);
t119 = sin(t123);
t120 = cos(t123);
t162 = -g(1) * t120 - g(2) * t119;
t130 = -pkin(2) - pkin(3);
t121 = qJDD(1) + qJDD(2);
t161 = t121 * pkin(2);
t160 = pkin(1) * qJD(1);
t124 = sin(qJ(4));
t122 = qJD(1) + qJD(2);
t125 = sin(qJ(2));
t146 = t125 * t160;
t94 = t122 * qJ(3) + t146;
t159 = t124 * t94;
t158 = pkin(1) * qJDD(1);
t117 = -qJDD(4) + t121;
t100 = t117 * MDP(10);
t157 = t121 * MDP(4) + t100;
t128 = cos(qJ(2));
t144 = qJD(2) * t160;
t156 = t125 * t158 + t128 * t144;
t155 = t120 * pkin(2) + t119 * qJ(3);
t154 = -t125 * t144 + t128 * t158;
t153 = qJD(2) * t125;
t152 = qJD(2) * t128;
t151 = qJD(4) * t124;
t127 = cos(qJ(4));
t150 = qJD(4) * t127;
t149 = -qJD(4) + t122;
t148 = t149 * t160;
t147 = pkin(1) * t153;
t145 = t128 * t160;
t111 = -t128 * pkin(1) - pkin(2);
t143 = t122 * t153;
t142 = -qJDD(3) + t154;
t141 = -t119 * pkin(2) + t120 * qJ(3);
t140 = t156 + t162;
t114 = t121 * qJ(3);
t116 = t122 * qJD(3);
t85 = t114 + t116 + t156;
t139 = qJD(3) - t145;
t138 = t127 * qJ(3) + t124 * t130;
t137 = t124 * qJ(3) - t127 * t130;
t136 = g(1) * t119 - g(2) * t120 + t154;
t135 = -qJDD(3) + t136;
t84 = t130 * t121 - t142;
t89 = t130 * t122 + t139;
t90 = -t119 * t124 - t120 * t127;
t91 = -t119 * t127 + t120 * t124;
t134 = g(1) * t90 + g(2) * t91 + t124 * t84 + t127 * t85 + t89 * t150;
t133 = g(1) * t91 - g(2) * t90 - t124 * t85 - t94 * t150 - t89 * t151;
t132 = t122 * t145 - t140;
t131 = t127 * t84 + t133;
t129 = cos(qJ(1));
t126 = sin(qJ(1));
t106 = t125 * pkin(1) + qJ(3);
t103 = -pkin(3) + t111;
t99 = pkin(1) * t152 + qJD(3);
t95 = t122 * t146;
t93 = -t122 * pkin(2) + t139;
t88 = -t142 - t161;
t1 = [qJDD(1) * MDP(1) + (g(1) * t126 - g(2) * t129) * MDP(2) + (g(1) * t129 + g(2) * t126) * MDP(3) + ((t121 * t128 - t143) * pkin(1) + t136) * MDP(5) + ((-t121 * t125 - t122 * t152) * pkin(1) - t140) * MDP(6) + (-pkin(1) * t143 + (pkin(2) - t111) * t121 + t135) * MDP(7) + (t106 * t121 + t99 * t122 + t162 + t85) * MDP(8) + (t85 * t106 + t94 * t99 + t88 * t111 + t93 * t147 - g(1) * (-t126 * pkin(1) + t141) - g(2) * (t129 * pkin(1) + t155)) * MDP(9) + ((-(-qJD(4) * t103 - t99) * t149 + t106 * t117) * t124 + (-(-qJD(4) * t106 + t147) * t149 - t103 * t117 - t84) * t127 - t133) * MDP(11) + ((t124 * t147 + t127 * t99) * t149 + (t124 * t103 + t127 * t106) * t117 + ((t127 * t103 - t124 * t106) * t149 - t159) * qJD(4) + t134) * MDP(12) + t157; (t136 + t95) * MDP(5) + t132 * MDP(6) + (t135 + t95 + 0.2e1 * t161) * MDP(7) + (0.2e1 * t114 + 0.2e1 * t116 - t132) * MDP(8) + (t85 * qJ(3) + t94 * qJD(3) - t88 * pkin(2) - g(1) * t141 - g(2) * t155 + (-t125 * t93 - t128 * t94) * t160) * MDP(9) + (-(-t124 * qJD(3) - qJD(4) * t138) * t149 + t137 * t117 + (-t124 * t128 + t125 * t127) * t148 - t131) * MDP(11) + (t127 * qJD(3) * t149 + t138 * t117 + (-t137 * t149 - t159) * qJD(4) - (t124 * t125 + t127 * t128) * t148 + t134) * MDP(12) + t157; -t121 * MDP(7) - t122 ^ 2 * MDP(8) + (-t94 * t122 - t135 - t161) * MDP(9) + (-t127 * MDP(11) + t124 * MDP(12)) * t117 - (MDP(11) * t124 + MDP(12) * t127) * t149 ^ 2; -t100 + ((-t124 * t89 - t127 * t94) * t149 + t131) * MDP(11) + (t94 * t151 - (t127 * t89 - t159) * t149 - t134) * MDP(12);];
tau  = t1;
