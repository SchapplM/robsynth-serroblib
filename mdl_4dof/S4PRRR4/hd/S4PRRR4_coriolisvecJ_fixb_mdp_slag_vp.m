% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR4
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:40
% EndTime: 2019-12-31 16:32:42
% DurationCPUTime: 0.44s
% Computational Cost: add. (262->80), mult. (687->133), div. (0->0), fcn. (426->4), ass. (0->53)
t102 = qJD(3) + qJD(4);
t138 = qJD(4) - t102;
t125 = (qJD(2) * qJD(3));
t143 = -2 * t125;
t106 = sin(qJ(3));
t142 = MDP(5) * t106;
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t108 = cos(qJ(3));
t89 = t105 * t108 + t107 * t106;
t141 = qJD(2) * t89;
t140 = (t106 ^ 2 - t108 ^ 2) * MDP(6);
t137 = pkin(5) + pkin(6);
t122 = qJD(2) * t137;
t83 = t108 * qJD(1) - t106 * t122;
t139 = MDP(10) * t106 + MDP(11) * t108;
t135 = qJD(3) * pkin(3);
t126 = t106 * qJD(1);
t94 = t137 * t108;
t84 = qJD(2) * t94 + t126;
t134 = t107 * t84;
t109 = qJD(3) ^ 2;
t133 = t106 * t109;
t132 = t108 * t109;
t128 = qJD(2) * t106;
t127 = qJD(2) * t108;
t124 = t106 * t135;
t123 = pkin(3) * t128;
t100 = -t108 * pkin(3) - pkin(2);
t121 = qJD(3) * t137;
t120 = t105 * t128;
t119 = t107 * t127;
t82 = t83 + t135;
t118 = -pkin(3) * t102 - t82;
t117 = t108 * t125;
t116 = pkin(2) * t143;
t88 = t105 * t106 - t107 * t108;
t80 = t83 * qJD(3);
t81 = (-t108 * t122 - t126) * qJD(3);
t87 = -t105 * t127 - t107 * t128;
t92 = t100 * qJD(2);
t113 = -t105 * t80 + t107 * t81 + t92 * t87;
t74 = qJD(4) * t119 - t102 * t120 + t107 * t117;
t85 = -t119 + t120;
t112 = -t87 * t85 * MDP(12) + t74 * MDP(14) + (-t85 ^ 2 + t87 ^ 2) * MDP(13) + (t85 * MDP(14) + (-t141 - t87) * MDP(15)) * t102;
t111 = t92 * t85 + (t138 * t84 - t81) * t105;
t77 = t102 * t89;
t93 = t137 * t106;
t91 = t108 * t121;
t90 = t106 * t121;
t76 = t102 * t88;
t75 = t77 * qJD(2);
t1 = [-t139 * t109 + (-t77 * MDP(17) + t76 * MDP(18)) * t102; 0.2e1 * t117 * t142 + t140 * t143 + MDP(7) * t132 - MDP(8) * t133 + (-pkin(5) * t132 + t106 * t116) * MDP(10) + (pkin(5) * t133 + t108 * t116) * MDP(11) + (t74 * t89 + t87 * t76) * MDP(12) + (-t74 * t88 - t89 * t75 + t76 * t85 + t87 * t77) * MDP(13) + (t100 * t75 + t92 * t77 + (qJD(2) * t88 + t85) * t124) * MDP(17) + (t100 * t74 - t92 * t76 + (-t87 + t141) * t124) * MDP(18) + (-t76 * MDP(14) - t77 * MDP(15) + (t105 * t90 - t107 * t91 + (t105 * t93 - t107 * t94) * qJD(4)) * MDP(17) + (t105 * t91 + t107 * t90 - (-t105 * t94 - t107 * t93) * qJD(4)) * MDP(18)) * t102; (-t85 * t123 - (-t105 * t83 - t134) * t102 + (t118 * t105 - t134) * qJD(4) + t113) * MDP(17) + (t87 * t123 + (t118 * qJD(4) + t83 * t102 - t80) * t107 + t111) * MDP(18) + t112 + (t139 * pkin(2) - t108 * t142 + t140) * qJD(2) ^ 2; (t113 + t138 * (-t105 * t82 - t134)) * MDP(17) + ((-t138 * t82 - t80) * t107 + t111) * MDP(18) + t112;];
tauc = t1;
