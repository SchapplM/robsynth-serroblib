% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:31
% EndTime: 2019-12-31 17:46:33
% DurationCPUTime: 0.57s
% Computational Cost: add. (321->86), mult. (699->140), div. (0->0), fcn. (418->6), ass. (0->50)
t111 = sin(pkin(8));
t112 = sin(pkin(7));
t113 = cos(pkin(8));
t152 = (t113 * MDP(10) - t111 * MDP(11) + MDP(7)) * t112;
t116 = sin(qJ(5));
t117 = cos(qJ(5));
t121 = t117 * t111 + t116 * t113;
t90 = t121 * qJD(1);
t139 = t111 ^ 2 + t113 ^ 2;
t114 = cos(pkin(7));
t101 = t114 * qJD(2) - qJD(4);
t97 = qJD(1) * t101;
t127 = t139 * t97;
t119 = qJD(1) ^ 2;
t150 = t139 * MDP(12) * t119;
t148 = qJ(2) * MDP(6) + t114 * MDP(8) + MDP(5);
t118 = -pkin(1) - pkin(2);
t140 = t114 * qJ(2) + t112 * t118;
t95 = -qJ(4) + t140;
t146 = pkin(6) - t95;
t137 = qJD(1) * qJ(2);
t99 = qJD(1) * t118 + qJD(2);
t87 = t112 * t99 + t114 * t137;
t142 = t116 * t111;
t141 = t117 * t113;
t136 = qJD(1) * t113;
t135 = qJD(2) * t112;
t134 = qJD(5) * t114;
t133 = t112 * qJD(5) ^ 2;
t130 = qJD(1) * t142;
t129 = t117 * t136;
t86 = -t112 * t137 + t114 * t99;
t126 = -t112 * qJ(2) + t114 * t118;
t124 = pkin(3) - t126;
t123 = t139 * (-qJD(1) * qJ(4) + t87);
t122 = t86 * t112 - t87 * t114;
t93 = -t141 + t142;
t80 = qJD(1) * pkin(3) + qJD(4) - t86;
t92 = t121 * qJD(5);
t120 = t123 + t135;
t98 = qJD(5) * t130;
t91 = t93 * qJD(5);
t89 = -t129 + t130;
t88 = t113 * pkin(4) + t124;
t85 = qJD(1) * t92;
t84 = -qJD(5) * t129 + t98;
t83 = t146 * t113;
t82 = t146 * t111;
t78 = pkin(4) * t136 + t80;
t1 = [(t95 * t127 + t123 * t101 + (qJD(1) * t124 + t80) * t135) * MDP(13) + (-t121 * t84 - t90 * t91) * MDP(14) + (-t121 * t85 + t84 * t93 + t91 * t89 - t90 * t92) * MDP(15) + t91 * qJD(5) * MDP(16) + t92 * qJD(5) * MDP(17) + (-t78 * t92 - t88 * t85 + ((-t116 * t82 + t117 * t83) * qJD(5) - t121 * t101) * qJD(5) + (-qJD(1) * t93 - t89) * t135) * MDP(19) + (t78 * t91 + t88 * t84 + ((-t116 * t83 - t117 * t82) * qJD(5) + t93 * t101) * qJD(5) - 0.2e1 * t90 * t135) * MDP(20) - 0.2e1 * MDP(12) * t127 + (((-t112 * t126 + t114 * t140) * qJD(1) - t122) * MDP(9) + (0.2e1 * t148 + 0.2e1 * t152) * qJD(1)) * qJD(2); t122 * qJD(1) * MDP(9) + t114 * t150 + (t112 * t127 + (-t80 * t112 - t114 * t120) * qJD(1)) * MDP(13) + (t114 * t85 + t93 * t133 + (t112 * t89 + t121 * t134) * qJD(1)) * MDP(19) + (-t114 * t84 + t121 * t133 + (t112 * t90 - t134 * t93) * qJD(1)) * MDP(20) + (-t148 - t152) * t119; (-t92 * MDP(19) + t91 * MDP(20)) * qJD(5); -t90 * qJD(5) * MDP(19) + (t89 * qJD(5) + t98) * MDP(20) - t150 + (t120 * MDP(13) + (-MDP(19) * t121 - MDP(20) * t141) * qJD(5)) * qJD(1); t90 * t89 * MDP(14) + (-t89 ^ 2 + t90 ^ 2) * MDP(15) + (t98 + (-t89 - t129) * qJD(5)) * MDP(16) + (-t121 * t97 + t78 * t90) * MDP(19) + (-t78 * t89 + t93 * t97) * MDP(20);];
tauc = t1;
