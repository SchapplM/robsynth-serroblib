% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:26
% EndTime: 2019-03-08 18:26:27
% DurationCPUTime: 0.42s
% Computational Cost: add. (183->101), mult. (643->168), div. (0->0), fcn. (394->4), ass. (0->60)
t100 = cos(pkin(4));
t135 = pkin(1) * t100;
t97 = sin(pkin(6));
t98 = sin(pkin(4));
t134 = t97 * t98;
t99 = cos(pkin(6));
t133 = t98 * t99;
t121 = t98 * qJD(1);
t110 = qJ(2) * t121;
t118 = qJD(1) * t100;
t112 = pkin(1) * t118;
t78 = t99 * t110 + t97 * t112;
t130 = qJ(2) * t98;
t132 = t99 * t130 + t97 * t135;
t93 = t97 ^ 2;
t95 = t99 ^ 2;
t131 = t93 + t95;
t126 = qJD(2) * t98;
t85 = t100 * qJD(3) + t99 * t126;
t129 = t100 * t85;
t82 = t85 * qJD(1);
t128 = t82 * t100;
t127 = MDP(12) * t98;
t94 = t98 ^ 2;
t125 = qJD(3) * t94;
t124 = t100 * qJ(3);
t123 = t97 * qJD(2);
t122 = t97 * qJD(3);
t111 = -pkin(1) * t99 - pkin(2);
t107 = -qJ(4) + t111;
t87 = t97 * t110;
t119 = qJD(3) + t87;
t70 = (pkin(3) * t134 + t107 * t100) * qJD(1) + t119;
t120 = qJD(3) + t70;
t117 = t100 * qJD(4);
t116 = -MDP(12) - MDP(8);
t115 = -MDP(14) + MDP(9);
t114 = MDP(10) + MDP(13);
t113 = 0.2e1 * qJD(2) * t94;
t109 = qJD(2) * t121;
t108 = -qJ(3) * t97 - pkin(1);
t106 = t111 * t100;
t86 = t97 * t109;
t81 = -qJD(1) * t117 + t86;
t105 = t81 * t97 + t82 * t99;
t80 = (-qJD(4) * t99 - t122) * t98;
t103 = (-pkin(2) * t99 + t108) * t98;
t102 = (-pkin(2) - qJ(4)) * t99 + t108;
t101 = qJD(1) ^ 2;
t96 = t100 ^ 2;
t90 = t97 * t130;
t84 = t98 * t123 - t117;
t77 = t99 * t112 - t87;
t76 = qJD(1) * t80;
t75 = -qJ(3) * t118 - t78;
t74 = qJD(1) * t103 + qJD(2);
t73 = qJD(1) * t106 + t119;
t72 = t102 * t121 + qJD(2);
t71 = qJD(4) + (pkin(3) * t133 + t124) * qJD(1) + t78;
t1 = [t131 * qJD(1) * MDP(6) * t113 + (-t77 * t97 + t78 * t99 + (t99 * t132 - t97 * (t99 * t135 - t90)) * qJD(1)) * MDP(7) * t126 + (t82 * t133 + (t93 * t113 + t85 * t133) * qJD(1)) * MDP(8) + 0.2e1 * (t100 * t126 - t99 * t125) * t97 * qJD(1) * MDP(9) + (t128 + (0.2e1 * t93 * t125 + t129) * qJD(1)) * MDP(10) + (-t82 * (-t124 - t132) - t75 * t85 + ((t73 * qJD(2) - t74 * qJD(3)) * t97 + ((t90 + t106) * t123 - t103 * t122) * qJD(1)) * t98) * MDP(11) + ((t84 * t97 + t85 * t99) * qJD(1) + t105) * t127 + (-t76 * t134 + t128 + (-t80 * t134 + t129) * qJD(1)) * MDP(13) + (-t76 * t133 - t81 * t100 + (-t100 * t84 - t80 * t133) * qJD(1)) * MDP(14) + (t72 * t80 + t81 * t90 + t70 * t84 + t82 * t132 + t71 * t85 + (t82 * qJ(3) + t81 * t107) * t100 + (t105 * pkin(3) + t76 * t102) * t98) * MDP(15) - 0.2e1 * (MDP(4) * t97 + MDP(5) * t99) * t100 * t109; -(MDP(6) - t116) * t131 * t94 * t101 + (((MDP(5) - t114) * t99 + (MDP(4) - t115) * t97) * t101 * t100 + ((-t78 * MDP(7) + t75 * MDP(11) + (-qJD(4) - t71) * MDP(15)) * t99 + (t77 * MDP(7) + (-qJD(3) - t73) * MDP(11) - t120 * MDP(15)) * t97) * qJD(1)) * t98; (MDP(11) + MDP(15)) * t86 + ((t100 * t75 + t74 * t134) * MDP(11) + (-t100 * t71 + t72 * t134 - t117) * MDP(15)) * qJD(1) + (t114 * (-t93 * t94 - t96) + (t116 * t98 * t100 + t115 * t97 * t94) * t99) * t101; ((-t94 * t95 - t96) * MDP(14) + (-t94 * t99 * MDP(13) + t100 * t127) * t97) * t101 + ((qJD(2) + t72) * t133 + t120 * t100) * MDP(15) * qJD(1);];
tauc  = t1;
