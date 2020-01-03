% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR5
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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:41
% DurationCPUTime: 0.25s
% Computational Cost: add. (216->61), mult. (377->95), div. (0->0), fcn. (152->4), ass. (0->41)
t79 = sin(qJ(3));
t105 = qJD(3) * t79;
t82 = -pkin(1) - pkin(2);
t70 = t82 * qJD(1) + qJD(2);
t81 = cos(qJ(3));
t98 = qJ(2) * qJD(1);
t92 = qJD(3) * t98;
t96 = qJD(1) * qJD(2);
t59 = t70 * t105 + t79 * t96 + t81 * t92;
t95 = qJD(1) - qJD(3);
t117 = -t59 - (t79 * t70 + t81 * t98) * t95;
t88 = t81 * qJ(2) + t79 * t82;
t116 = t59 + (t79 * qJD(2) + t88 * qJD(3)) * t95;
t80 = cos(qJ(4));
t101 = t80 * MDP(16);
t78 = sin(qJ(4));
t115 = t78 * MDP(15) + t101;
t114 = (t78 ^ 2 - t80 ^ 2) * MDP(11);
t113 = 0.2e1 * qJD(2);
t112 = t95 * pkin(3);
t109 = t81 * t70;
t108 = qJD(3) * t109 + t81 * t96;
t106 = MDP(16) * t78;
t104 = t95 * qJD(4) * t114;
t102 = t80 * MDP(10);
t83 = qJD(4) ^ 2;
t100 = t83 * MDP(12);
t99 = t83 * MDP(13);
t97 = MDP(16) * qJD(4);
t94 = t95 * t102;
t63 = -t79 * t98 + t109;
t60 = -t63 + t112;
t91 = t60 + t63 + t112;
t90 = t95 * MDP(15);
t87 = -t79 * qJ(2) + t81 * t82;
t61 = t81 * qJD(2) + t87 * qJD(3);
t89 = -(pkin(3) - t87) * t95 - t60 - t61;
t86 = pkin(6) * t83 - t117;
t85 = (-pkin(6) + t88) * t83 - t116;
t58 = -t79 * t92 + t108;
t1 = [t116 * MDP(8) + (t61 * t95 + t108) * MDP(9) - 0.2e1 * t104 + (MDP(5) * t113 + (MDP(6) * t113 - MDP(9) * t105) * qJ(2)) * qJD(1) + (-t85 * MDP(15) + t89 * t97 - t100) * t80 + (t99 + t85 * MDP(16) + (t89 * MDP(15) + 0.2e1 * t94) * qJD(4)) * t78; (-qJ(2) * MDP(6) - MDP(5)) * qJD(1) ^ 2 + (-MDP(15) * t80 + t106) * t83 * t79 + (t95 * t101 + t78 * t90) * t81 * qJD(4) - ((t95 * MDP(9) - qJD(4) * t115) * t81 + (t80 * t90 + (MDP(8) - t106) * t95) * t79) * t95; t117 * MDP(8) + (-t63 * t95 - t58) * MDP(9) + 0.2e1 * t104 + (-t86 * MDP(15) + t91 * t97 + t100) * t80 + (-t99 + t86 * MDP(16) + (t91 * MDP(15) - 0.2e1 * t94) * qJD(4)) * t78; (-t78 * t102 + t114) * t95 ^ 2 + t115 * (t60 * t95 - t58);];
tauc = t1;
