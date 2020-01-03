% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:02
% EndTime: 2019-12-31 16:45:04
% DurationCPUTime: 0.50s
% Computational Cost: add. (470->122), mult. (1270->165), div. (0->0), fcn. (842->4), ass. (0->57)
t124 = sin(pkin(6));
t125 = cos(pkin(6));
t152 = (t124 ^ 2 + t125 ^ 2) * (qJ(2) * MDP(7) + MDP(6));
t126 = sin(qJ(3));
t127 = cos(qJ(3));
t108 = t127 * t124 + t126 * t125;
t103 = t108 * qJD(1);
t139 = MDP(13) + MDP(15);
t150 = t103 ^ 2;
t149 = pkin(5) + qJ(2);
t113 = t149 * t125;
t110 = qJD(1) * t113;
t147 = t126 * t110;
t112 = t149 * t124;
t109 = qJD(1) * t112;
t91 = -t127 * t109 - t147;
t146 = qJD(4) - t91;
t144 = qJD(1) * t126;
t143 = qJD(1) * t127;
t142 = qJD(3) * t126;
t141 = qJD(3) * t127;
t140 = qJD(1) * qJD(2);
t138 = -MDP(14) + MDP(17);
t121 = -t125 * pkin(2) - pkin(1);
t137 = t124 * t144;
t136 = t125 * t143;
t135 = t126 * t140;
t134 = t127 * t140;
t133 = t91 + t147;
t80 = -t109 * t142 + t110 * t141 + t124 * t134 + t125 * t135;
t117 = qJD(3) * t136;
t95 = qJD(3) * t137 - t117;
t106 = t108 * qJD(3);
t96 = qJD(1) * t106;
t132 = t96 * pkin(3) + t95 * qJ(4);
t92 = -t126 * t109 + t127 * t110;
t131 = -t127 * t112 - t126 * t113;
t94 = -t126 * t112 + t127 * t113;
t107 = t126 * t124 - t127 * t125;
t101 = -t136 + t137;
t111 = t121 * qJD(1) + qJD(2);
t82 = t101 * pkin(3) - t103 * qJ(4) + t111;
t130 = -t82 * t103 - t80;
t129 = t109 * t141 + t124 * t135 - t125 * t134;
t105 = t124 * t142 - t125 * t141;
t100 = t101 ^ 2;
t90 = t107 * pkin(3) - t108 * qJ(4) + t121;
t89 = t103 * pkin(3) + t101 * qJ(4);
t88 = qJD(3) * qJ(4) + t92;
t87 = t117 + (t101 - t137) * qJD(3);
t85 = -qJD(3) * pkin(3) + t146;
t84 = t108 * qJD(2) + t94 * qJD(3);
t83 = -t107 * qJD(2) + t131 * qJD(3);
t81 = t106 * pkin(3) + t105 * qJ(4) - t108 * qJD(4);
t79 = (qJD(4) - t147) * qJD(3) - t129;
t78 = -t103 * qJD(4) + t132;
t1 = [(-t103 * t105 - t95 * t108) * MDP(8) + (t105 * t101 - t103 * t106 + t95 * t107 - t108 * t96) * MDP(9) + (t111 * t106 + t121 * t96) * MDP(13) + (-t111 * t105 - t121 * t95) * MDP(14) + (t81 * t101 + t82 * t106 + t78 * t107 + t90 * t96) * MDP(15) + (-t83 * t101 + t84 * t103 - t85 * t105 - t88 * t106 - t79 * t107 + t80 * t108 + t131 * t95 - t94 * t96) * MDP(16) + (-t81 * t103 + t82 * t105 - t78 * t108 + t90 * t95) * MDP(17) + (-t131 * t80 + t78 * t90 + t79 * t94 + t82 * t81 + t88 * t83 + t85 * t84) * MDP(18) + (-t105 * MDP(10) - t106 * MDP(11) + t138 * t83 - t139 * t84) * qJD(3) + 0.2e1 * t140 * t152; (-t100 - t150) * MDP(16) + (t88 * t101 + (-qJD(4) - t85) * t103 + t132) * MDP(18) + t138 * (-t117 + (t101 + t137) * qJD(3)) - qJD(1) ^ 2 * t152 + 0.2e1 * t139 * t103 * qJD(3); (-t100 + t150) * MDP(9) + t87 * MDP(10) + (-t111 * t103 - t80) * MDP(13) + t129 * MDP(14) + t130 * MDP(15) + (pkin(3) * t95 - t96 * qJ(4) - (-t88 + t92) * t103) * MDP(16) + (t89 * t103 - t129) * MDP(17) + (-t80 * pkin(3) + t79 * qJ(4) + t146 * t88 - t82 * t89 - t85 * t92) * MDP(18) + (t103 * MDP(8) + t111 * MDP(14) - t89 * MDP(15) + (t85 - t146) * MDP(16) - t82 * MDP(17)) * t101 + ((-t124 * t143 - t125 * t144 + t103) * MDP(11) + t133 * MDP(14) + (0.2e1 * qJD(4) - t133) * MDP(17) + t139 * t92) * qJD(3); t103 * t101 * MDP(15) + t87 * MDP(16) + (-qJD(3) ^ 2 - t150) * MDP(17) + (-t88 * qJD(3) - t130) * MDP(18);];
tauc = t1;
