% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:25
% EndTime: 2019-12-31 17:32:26
% DurationCPUTime: 0.41s
% Computational Cost: add. (208->72), mult. (557->118), div. (0->0), fcn. (374->6), ass. (0->41)
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t105 = t82 ^ 2 + t83 ^ 2;
t88 = qJD(3) ^ 2;
t110 = t105 * MDP(8) * t88;
t87 = cos(qJ(3));
t99 = t87 * qJD(2);
t93 = qJD(4) - t99;
t86 = cos(qJ(5));
t107 = t86 * t83;
t84 = sin(qJ(5));
t109 = t84 * t82;
t68 = -t107 + t109;
t106 = pkin(6) + qJ(4);
t103 = qJD(3) * pkin(3);
t69 = t82 * t86 + t83 * t84;
t64 = qJD(3) * t69;
t102 = qJD(5) * t87;
t85 = sin(qJ(3));
t101 = t85 * qJD(2);
t100 = t85 * qJD(5) ^ 2;
t97 = qJD(3) * t109;
t96 = qJD(3) * t107;
t79 = -pkin(4) * t83 - pkin(3);
t71 = (qJD(4) + t99) * qJD(3);
t95 = t105 * t71;
t92 = t105 * (qJD(3) * qJ(4) + t101);
t66 = t69 * qJD(5);
t91 = t69 * t102;
t90 = t68 * t102;
t89 = -t92 + t101;
t75 = qJD(5) * t96;
t74 = t93 - t103;
t73 = t106 * t83;
t72 = t106 * t82;
t67 = qJD(3) * t79 + t93;
t65 = t68 * qJD(5);
t62 = -t96 + t97;
t59 = qJD(3) * t66;
t58 = -qJD(5) * t97 + t75;
t1 = [(MDP(15) * t66 - MDP(16) * t65) * qJD(5); t87 * t110 + (t85 * t95 + (t74 * t85 - t87 * t89) * qJD(3)) * MDP(9) + (-t87 * t59 + t68 * t100 + (t85 * t62 - t91) * qJD(3)) * MDP(15) + (-t87 * t58 + t69 * t100 + (t85 * t64 + t90) * qJD(3)) * MDP(16) + (-t87 * MDP(5) + (-MDP(6) * t83 + MDP(7) * t82 - MDP(4)) * t85) * t88; (qJD(3) * t105 * t93 + t95) * MDP(8) + (t92 * qJD(4) + qJ(4) * t95 + (-t92 * t87 + (-t74 - t103) * t85) * qJD(2)) * MDP(9) + (t58 * t69 - t64 * t65) * MDP(10) + (-t58 * t68 - t59 * t69 + t62 * t65 - t64 * t66) * MDP(11) - t65 * qJD(5) * MDP(12) - t66 * qJD(5) * MDP(13) + (t79 * t59 + t67 * t66 + ((t72 * t84 - t73 * t86) * qJD(5) - t69 * qJD(4)) * qJD(5) + (t91 + (qJD(3) * t68 - t62) * t85) * qJD(2)) * MDP(15) + (t79 * t58 - t67 * t65 + ((t72 * t86 + t73 * t84) * qJD(5) + t68 * qJD(4)) * qJD(5) - t90 * qJD(2)) * MDP(16); t64 * qJD(5) * MDP(15) + (-qJD(5) * t62 + t75) * MDP(16) - t110 + (t89 * MDP(9) + (MDP(15) * t69 - MDP(16) * t109) * qJD(5)) * qJD(3); t64 * t62 * MDP(10) + (-t62 ^ 2 + t64 ^ 2) * MDP(11) + (t75 + (t62 - t97) * qJD(5)) * MDP(12) + (-t67 * t64 - t69 * t71) * MDP(15) + (t67 * t62 + t68 * t71) * MDP(16);];
tauc = t1;
