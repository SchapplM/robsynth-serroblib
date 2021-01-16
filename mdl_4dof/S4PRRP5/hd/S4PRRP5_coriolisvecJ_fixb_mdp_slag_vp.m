% Calculate Coriolis joint torque vector for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:10
% EndTime: 2021-01-14 22:36:12
% DurationCPUTime: 0.47s
% Computational Cost: add. (237->99), mult. (611->144), div. (0->0), fcn. (301->4), ass. (0->51)
t90 = sin(qJ(2));
t93 = qJD(3) ^ 2;
t94 = qJD(2) ^ 2;
t126 = (t93 + t94) * t90;
t89 = sin(qJ(3));
t87 = t89 ^ 2;
t91 = cos(qJ(3));
t88 = t91 ^ 2;
t125 = (t87 - t88) * MDP(6);
t101 = (t87 + t88) * MDP(14);
t124 = t91 * pkin(3);
t123 = qJ(4) + pkin(5);
t116 = qJD(3) * pkin(3);
t110 = t90 * qJD(1);
t82 = qJD(2) * pkin(5) + t110;
t99 = qJ(4) * qJD(2) + t82;
t72 = t99 * t89;
t71 = -t72 + t116;
t122 = t71 + t72;
t106 = qJD(2) * qJD(3);
t77 = t89 * pkin(3) * t106 + qJD(2) * t110;
t118 = pkin(3) * MDP(15);
t117 = qJD(2) * pkin(2);
t115 = t91 * MDP(5);
t86 = -pkin(2) - t124;
t114 = qJD(2) * t86;
t113 = qJD(3) * t89;
t112 = qJD(3) * t91;
t111 = t87 * MDP(13);
t109 = t91 * MDP(12);
t92 = cos(qJ(2));
t108 = t92 * qJD(1);
t107 = qJ(4) * qJD(3);
t105 = MDP(10) + MDP(12);
t104 = MDP(11) + MDP(13);
t102 = qJD(3) * t123;
t83 = -t108 - t117;
t100 = -t83 - t108;
t98 = -0.2e1 * t92 * t106;
t97 = qJD(4) + t108;
t73 = t99 * t91;
t96 = t71 * t89 - t73 * t91;
t76 = qJD(4) - t108 + t114;
t95 = -t76 - t97;
t79 = t123 * t91;
t78 = t123 * t89;
t75 = -t89 * qJD(4) - t91 * t102;
t74 = t91 * qJD(4) - t89 * t102;
t68 = -t82 * t112 + (-t91 * t107 - t97 * t89) * qJD(2);
t67 = -t82 * t113 + (-t89 * t107 + t97 * t91) * qJD(2);
t1 = [t105 * (-t91 * t126 + t89 * t98) + t104 * (t89 * t126 + t91 * t98) + ((-t96 * qJD(2) - t77) * MDP(15) + (-MDP(4) + t101) * t94) * t92 + (-t94 * MDP(3) + (qJD(2) * t76 - t71 * t112 - t73 * t113 + t67 * t91 - t68 * t89) * MDP(15)) * t90; (t67 * t79 - t68 * t78 + t71 * t75 + t73 * t74 + t77 * t86) * MDP(15) + (-t77 * MDP(12) + (qJD(2) * t74 + t67) * MDP(14) + (-pkin(5) * MDP(10) + MDP(7)) * t93) * t91 + (t77 * MDP(13) + (-qJD(2) * t75 - t68) * MDP(14) + (pkin(5) * MDP(11) - MDP(8)) * t93) * t89 + ((-t76 * t90 + t96 * t92) * MDP(15) + (-t92 * t101 + (-t89 * MDP(13) + t109) * t90) * qJD(2)) * qJD(1) + (t75 * MDP(12) - t74 * MDP(13) + (pkin(3) * t111 - 0.2e1 * t125) * qJD(2) + ((t83 - t117) * MDP(11) + (t76 + t114) * MDP(13) + (qJD(2) * t78 - t71) * MDP(14)) * t91 + (t83 * MDP(10) - t73 * MDP(14) + (MDP(12) + t118) * t76 + (0.2e1 * t115 - pkin(2) * MDP(10) + (t86 - t124) * MDP(12) - t79 * MDP(14)) * qJD(2)) * t89 + (t104 * t91 + t105 * t89) * t108) * qJD(3); (t68 * pkin(3) + t122 * t73) * MDP(15) + ((-t91 * t82 + t73) * MDP(12) + (t89 * t82 - t72) * MDP(13)) * qJD(3) + (-t89 * t115 + t125 + (t89 * t109 - t111) * pkin(3)) * t94 + ((t100 * MDP(10) + t95 * MDP(12) + MDP(13) * t107 - t76 * t118) * t89 + (t100 * MDP(11) - MDP(12) * t107 + t95 * MDP(13) + (-t116 + t122) * MDP(14)) * t91) * qJD(2); t77 * MDP(15) - t94 * t101 + (t96 * MDP(15) + 0.2e1 * (t89 * MDP(12) + t91 * MDP(13)) * qJD(3)) * qJD(2);];
tauc = t1;
