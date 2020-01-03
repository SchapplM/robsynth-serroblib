% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:37
% EndTime: 2019-12-31 19:26:38
% DurationCPUTime: 0.36s
% Computational Cost: add. (319->94), mult. (691->137), div. (0->0), fcn. (344->6), ass. (0->56)
t102 = qJD(1) + qJD(2);
t105 = sin(pkin(8));
t106 = cos(pkin(8));
t108 = sin(qJ(2));
t139 = pkin(1) * qJD(1);
t124 = t108 * t139;
t110 = cos(qJ(2));
t123 = t110 * t139;
t90 = t102 * pkin(2) + t123;
t79 = t105 * t90 + t106 * t124;
t74 = t102 * qJ(4) + t79;
t138 = pkin(1) * qJD(2);
t122 = qJD(1) * t138;
t116 = t110 * t122;
t117 = t108 * t122;
t81 = t105 * t116 + t106 * t117;
t120 = t74 * t102 - t81;
t107 = sin(qJ(5));
t109 = cos(qJ(5));
t142 = (t107 ^ 2 - t109 ^ 2) * MDP(12);
t101 = t102 ^ 2;
t128 = t109 * qJD(5);
t82 = -t105 * t117 + t106 * t116;
t77 = t102 * qJD(4) + t82;
t141 = t77 * t107 + t74 * t128;
t140 = 0.2e1 * (-MDP(11) * t107 * t128 + qJD(5) * t142) * t102;
t111 = qJD(5) ^ 2;
t121 = -t106 * pkin(2) - pkin(3);
t137 = t111 * (-pkin(7) + t121);
t100 = t110 * pkin(1) + pkin(2);
t135 = t105 * t100;
t134 = t106 * t108;
t133 = t108 * MDP(5);
t95 = t105 * t124;
t88 = t106 * t123 - t95;
t132 = qJD(4) - t88;
t129 = t109 * MDP(16);
t127 = t111 * MDP(13);
t126 = t111 * MDP(14);
t125 = -qJD(1) - t102;
t97 = t105 * t108 * pkin(1);
t85 = pkin(1) * t134 + qJ(4) + t135;
t112 = pkin(1) * (t105 * t110 + t134);
t87 = qJD(2) * t112;
t119 = t102 * t85 + t87;
t78 = t106 * t90 - t95;
t118 = t106 * t100 - t97;
t115 = -pkin(3) - t118;
t96 = qJD(2) * t97;
t83 = t106 * t110 * t138 + qJD(4) - t96;
t113 = t83 * t102 - t111 * (-pkin(7) + t115);
t99 = t105 * pkin(2) + qJ(4);
t86 = qJD(1) * t112;
t76 = t77 * t109;
t73 = -t102 * pkin(3) + qJD(4) - t78;
t1 = [(-t81 * t118 + t82 * t135 - t78 * t87 - t79 * t96) * MDP(7) + (t87 * t102 + t81) * MDP(8) + ((qJD(4) + t83) * t102 + t82) * MDP(9) + (t81 * t115 + t73 * t87 + t74 * t83 + t77 * t85) * MDP(10) + t141 * MDP(16) + t76 * MDP(17) + (t119 * MDP(16) * qJD(5) + t113 * MDP(17) - t126) * t109 + (t82 * MDP(7) * t134 + (t125 * t133 + (t106 * t79 * MDP(7) + t125 * MDP(6)) * t110) * qJD(2)) * pkin(1) + (-t127 + t113 * MDP(16) + (-t119 - t74) * MDP(17) * qJD(5)) * t107 + t140; (t78 * t86 - t79 * t88 + (t82 * t105 - t81 * t106) * pkin(2)) * MDP(7) + (-t86 * t102 + t81) * MDP(8) + ((0.2e1 * qJD(4) - t88) * t102 + t82) * MDP(9) + (t81 * t121 + t132 * t74 - t73 * t86 + t77 * t99) * MDP(10) - t107 * t127 - t109 * t126 + (-t86 * t128 - t107 * t137 + (t132 * t107 + t99 * t128) * t102 + t141) * MDP(16) + (t76 + (t132 * t102 - t137) * t109 + (-t102 * t99 - t74 + t86) * t107 * qJD(5)) * MDP(17) + t140 + (t110 * MDP(6) + t133) * (-qJD(2) + t102) * t139; (t107 * MDP(17) - t129) * t111; -t101 * MDP(9) - t120 * MDP(10) + (t107 * MDP(16) + t109 * MDP(17)) * (-t101 - t111); -t120 * t129 - t101 * t142 + (t109 * t101 * MDP(11) + t120 * MDP(17)) * t107;];
tauc = t1;
