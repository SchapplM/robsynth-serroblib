% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:15
% EndTime: 2019-12-05 15:07:17
% DurationCPUTime: 0.47s
% Computational Cost: add. (338->79), mult. (878->128), div. (0->0), fcn. (586->6), ass. (0->53)
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t142 = -t100 * t97 + t102 * t98;
t82 = t142 * qJD(1);
t109 = t100 * t98 + t102 * t97;
t85 = t109 * qJD(3);
t135 = -qJ(5) - pkin(6);
t79 = qJD(1) * t85;
t83 = t109 * qJD(1);
t141 = t83 * qJD(3) - t79;
t101 = cos(qJ(4));
t112 = -t135 * qJD(3) + t83;
t99 = sin(qJ(4));
t72 = t101 * qJD(2) - t112 * t99;
t73 = t99 * qJD(2) + t112 * t101;
t140 = qJD(4) * t73;
t95 = t99 ^ 2;
t96 = t101 ^ 2;
t139 = (t95 - t96) * MDP(7);
t138 = t99 * MDP(11) + t101 * MDP(12);
t132 = t95 + t96;
t137 = t132 * MDP(13);
t103 = qJD(4) ^ 2;
t84 = t142 * qJD(3);
t78 = qJD(1) * t84;
t111 = qJD(3) * qJD(5) + t78;
t70 = -t111 * t99 - t140;
t136 = (t70 + t140) * MDP(14) - t103 * MDP(12);
t130 = qJD(4) * pkin(4);
t71 = t72 + t130;
t134 = t71 - t72;
t131 = qJD(3) * pkin(3);
t122 = t99 * qJD(4);
t118 = qJD(3) * qJD(4);
t116 = t99 * t118;
t75 = pkin(4) * t116 + t79;
t117 = -t101 * pkin(4) - pkin(3);
t115 = qJD(4) * t135;
t110 = t101 * t73 - t71 * t99;
t107 = pkin(6) * t103 - t141;
t76 = -t82 - t131;
t106 = qJD(4) * (t76 + t82 - t131);
t69 = t72 * qJD(4) + t111 * t101;
t105 = -t103 * MDP(11) + (-qJD(4) * t71 + t69) * MDP(14);
t104 = qJD(3) ^ 2;
t89 = t135 * t101;
t88 = t135 * t99;
t81 = -t99 * qJD(5) + t101 * t115;
t80 = t101 * qJD(5) + t99 * t115;
t74 = t117 * qJD(3) + qJD(5) - t82;
t1 = [(-t142 * t75 + t74 * t85) * MDP(14) + (t110 * MDP(14) - t138 * qJD(4)) * t84 + (t105 * t101 - t136 * t99) * t109 + (-t85 * MDP(4) + (-t101 * t85 - t122 * t142) * MDP(11) + (-qJD(4) * t101 * t142 + t85 * t99) * MDP(12) + (-MDP(5) + t137) * t84) * qJD(3); t136 * t101 + t105 * t99; t141 * MDP(4) + 0.2e1 * t101 * MDP(6) * t116 - 0.2e1 * t118 * t139 + (-t107 * t101 + t99 * t106) * MDP(11) + (t101 * t106 + t107 * t99) * MDP(12) + (t69 * t101 - t70 * t99 + (-t101 * t71 - t73 * t99) * qJD(4) + (t101 * t80 - t81 * t99 - t132 * t82 + (-t101 * t88 + t89 * t99) * qJD(4)) * qJD(3)) * MDP(13) + (-t69 * t89 + t70 * t88 + t75 * t117 + (pkin(4) * t122 - t83) * t74 + (-t101 * t82 + t80) * t73 + (t99 * t82 + t81) * t71) * MDP(14) + (t101 * MDP(8) - t99 * MDP(9)) * t103; t104 * t139 + (t134 * t73 + (-qJD(3) * t74 * t99 + t70) * pkin(4)) * MDP(14) + t138 * (-t76 * qJD(3) - t78) + (-t99 * t104 * MDP(6) + (-t130 + t134) * qJD(3) * MDP(13)) * t101; (-t110 * qJD(3) + t75) * MDP(14) - t104 * t137;];
tauc = t1;
