% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:47
% EndTime: 2019-12-31 16:33:48
% DurationCPUTime: 0.27s
% Computational Cost: add. (204->67), mult. (457->112), div. (0->0), fcn. (284->6), ass. (0->45)
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t138 = t96 * t99 * MDP(8) - (t96 ^ 2 - t99 ^ 2) * MDP(9);
t93 = qJD(2) + qJD(3);
t97 = sin(qJ(3));
t126 = qJD(3) * t97;
t100 = cos(qJ(3));
t101 = cos(qJ(2));
t98 = sin(qJ(2));
t85 = t100 * t98 + t97 * t101;
t137 = -pkin(2) * t126 + t85 * qJD(1);
t135 = t96 * MDP(13) + t99 * MDP(14);
t134 = t93 * pkin(3);
t127 = qJD(1) * t98;
t119 = t101 * qJD(1);
t88 = qJD(2) * pkin(2) + t119;
t133 = (t100 * t127 + t97 * t88) * t93;
t124 = qJD(4) * t99;
t105 = t85 * qJD(2);
t120 = qJD(3) * t100;
t112 = t98 * t120;
t114 = t88 * t126;
t72 = t114 + (t105 + t112) * qJD(1);
t113 = t97 * t127;
t79 = t100 * t88 - t113;
t77 = -t79 - t134;
t132 = t77 * t124 + t72 * t96;
t131 = t93 * t113;
t102 = qJD(4) ^ 2;
t129 = t102 * t96;
t128 = t102 * t99;
t125 = qJD(4) * t96;
t123 = t100 * t101;
t118 = t102 * MDP(11);
t109 = pkin(6) * t102 - t133;
t108 = qJD(4) * (t79 - t134);
t107 = -t97 * t98 + t123;
t106 = -t88 * t120 + t131;
t104 = MDP(10) * t128 + (-MDP(6) * t112 + (-t85 * MDP(6) - MDP(7) * t123) * qJD(2)) * qJD(1) + 0.2e1 * t138 * qJD(4) * t93;
t91 = -t100 * pkin(2) - pkin(3);
t90 = t97 * pkin(2) + pkin(6);
t75 = t77 * t125;
t74 = t85 * qJD(3) + t105;
t73 = t93 * t107;
t1 = [(-t73 * t125 - t85 * t128) * MDP(13) + (-t73 * t124 + t85 * t129) * MDP(14) + (-t98 * MDP(3) - t101 * MDP(4)) * qJD(2) ^ 2 + (-t74 * MDP(6) - t73 * MDP(7) + (-t107 * t125 - t74 * t99) * MDP(13) + (-t107 * t124 + t74 * t96) * MDP(14)) * t93; -MDP(6) * t114 + t106 * MDP(7) - t96 * t118 + (-t90 * t128 - t72 * t99 + t75) * MDP(13) + (t90 * t129 + t132) * MDP(14) + (t137 * MDP(6) + (t91 * t125 + t137 * t99) * MDP(13) + (t91 * t124 - t137 * t96) * MDP(14)) * t93 + t104 + (t93 * MDP(7) + t135 * qJD(4)) * (-pkin(2) * t120 + t107 * qJD(1)); (-t114 + t133) * MDP(6) + (t79 * t93 + t106) * MDP(7) + t75 * MDP(13) + t132 * MDP(14) + ((-t109 - t72) * MDP(13) + MDP(14) * t108) * t99 + (MDP(13) * t108 + t109 * MDP(14) - t118) * t96 + t104; -t138 * t93 ^ 2 + t135 * (-t77 * t93 - (qJD(2) * t119 + qJD(3) * t88) * t100 + t131);];
tauc = t1;
