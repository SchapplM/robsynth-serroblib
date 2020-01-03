% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:30
% EndTime: 2019-12-31 16:56:32
% DurationCPUTime: 0.85s
% Computational Cost: add. (376->141), mult. (855->232), div. (0->0), fcn. (468->4), ass. (0->68)
t120 = 2 * qJD(1);
t101 = cos(qJ(3));
t97 = t101 ^ 2;
t99 = sin(qJ(3));
t148 = (t99 ^ 2 - t97) * MDP(8);
t147 = (MDP(6) * qJ(2) + MDP(5));
t102 = (-pkin(1) - pkin(5));
t98 = sin(qJ(4));
t131 = qJD(4) * t98;
t114 = t101 * t131;
t100 = cos(qJ(4));
t122 = t100 * qJD(3);
t94 = qJD(4) * t122;
t78 = t94 + (-t99 * t122 - t114) * qJD(1);
t146 = t78 * t98;
t92 = t102 * qJD(1) + qJD(2);
t84 = -qJD(3) * pkin(3) - t101 * t92;
t145 = t84 * t98;
t124 = qJD(1) * t101;
t117 = t98 * t124;
t86 = t117 - t122;
t134 = qJD(1) * t99;
t93 = qJD(4) + t134;
t144 = t86 * t93;
t133 = qJD(3) * t98;
t88 = t100 * t124 + t133;
t143 = t88 * t93;
t142 = t93 * t98;
t141 = t98 * t99;
t139 = t100 * t84;
t138 = t100 * t93;
t137 = t101 * t98;
t136 = t102 * t98;
t132 = qJD(3) * t99;
t130 = t100 * t101;
t129 = t100 * t102;
t128 = t101 * t102;
t104 = qJD(1) ^ 2;
t127 = t101 * t104;
t103 = qJD(3) ^ 2;
t126 = t102 * t103;
t125 = -t103 - t104;
t123 = qJD(3) * t101;
t121 = qJD(1) * qJD(3);
t119 = t93 * t131;
t116 = t92 * t123;
t115 = qJD(4) * t138;
t113 = t99 * t121;
t112 = MDP(18) * t123;
t111 = t93 + t134;
t110 = -t88 + t133;
t109 = t86 + t122;
t108 = pkin(3) * t101 + pkin(6) * t99;
t90 = pkin(3) * t99 - pkin(6) * t101 + qJ(2);
t82 = t90 * qJD(1);
t83 = qJD(3) * pkin(6) + t92 * t99;
t77 = t100 * t83 + t82 * t98;
t76 = t100 * t82 - t83 * t98;
t107 = qJD(1) * t97 - t93 * t99;
t106 = t100 * t90 - t99 * t136;
t105 = t99 * t129 + t90 * t98;
t85 = t108 * qJD(3) + qJD(2);
t91 = t98 * t113;
t89 = t108 * qJD(1);
t81 = t85 * qJD(1);
t80 = t100 * t81;
t79 = qJD(4) * t88 - t91;
t1 = [-0.2e1 * t101 * MDP(7) * t113 + 0.2e1 * t121 * t148 + (-t99 * t126 + (qJ(2) * t123 + qJD(2) * t99) * t120) * MDP(12) + (-t101 * t126 + (-qJ(2) * t132 + qJD(2) * t101) * t120) * MDP(13) + (-t88 * t114 + (t101 * t78 - t88 * t132) * t100) * MDP(14) + ((t100 * t86 + t88 * t98) * t132 + (-t100 * t79 - t146 + (-t100 * t88 + t86 * t98) * qJD(4)) * t101) * MDP(15) + (-t93 * t114 + t78 * t99 + (t107 * t100 + t101 * t88) * qJD(3)) * MDP(16) + (-t101 * t115 - t79 * t99 + (-t101 * t86 - t107 * t98) * qJD(3)) * MDP(17) + t111 * t112 + (t85 * t138 - t79 * t128 + t80 * t99 + (-t105 * t93 + t84 * t130 - t77 * t99) * qJD(4) + ((t102 * t86 - t145) * t99 + (t106 * qJD(1) - t93 * t136 + t76) * t101) * qJD(3)) * MDP(19) + (-t78 * t128 + (-t81 * t99 - t85 * t93) * t98 + (-t106 * t93 - t84 * t137 - t76 * t99) * qJD(4) + ((t102 * t88 - t139) * t99 + (-t105 * qJD(1) - t93 * t129 - t77) * t101) * qJD(3)) * MDP(20) + (t147 * qJD(2) * t120) + (-t101 * MDP(10) - t99 * MDP(9)) * t103; -(t147 * t104) + (-t100 * MDP(19) + t98 * MDP(20)) * t93 * qJD(1) + (t125 * MDP(12) + (qJD(3) * t86 - t115) * MDP(19) + (qJD(3) * t88 + t119) * MDP(20)) * t99 + (t125 * MDP(13) - t79 * MDP(19) - t78 * MDP(20) + (-MDP(19) * t98 - MDP(20) * t100) * qJD(3) * t111) * t101; t99 * MDP(7) * t127 - t104 * t148 + (t88 * t138 + t146) * MDP(14) + ((-t79 - t143) * t98 + (t78 - t144) * t100) * MDP(15) + (t115 + (t110 * t101 + t99 * t138) * qJD(1)) * MDP(16) + (-t119 + (t109 * t101 - t93 * t141) * qJD(1)) * MDP(17) - t93 * MDP(18) * t124 + (-t89 * t138 - pkin(3) * t79 + (-t109 * t99 + t93 * t137) * t92 + (-pkin(6) * t138 + t145) * qJD(4) + (t84 * t141 + (-pkin(6) * t133 - t76) * t101) * qJD(1)) * MDP(19) + (t89 * t142 - pkin(3) * t78 + (t110 * t99 + t93 * t130) * t92 + (pkin(6) * t142 + t139) * qJD(4) + (t77 * t101 + (-pkin(6) * t123 + t84 * t99) * t100) * qJD(1)) * MDP(20) + (t104 * t99 * MDP(13) - MDP(12) * t127) * qJ(2); t88 * t86 * MDP(14) + (-t86 ^ 2 + t88 ^ 2) * MDP(15) + (-t100 * t113 + t144 + t94) * MDP(16) + (t91 + t143) * MDP(17) + qJD(1) * t112 + (-t98 * t116 + t77 * t93 - t84 * t88 + t80) * MDP(19) + (-t100 * t116 + t76 * t93 - t81 * t98 + t84 * t86) * MDP(20) + (-MDP(16) * t117 - t88 * MDP(17) - t77 * MDP(19) - t76 * MDP(20)) * qJD(4);];
tauc = t1;
