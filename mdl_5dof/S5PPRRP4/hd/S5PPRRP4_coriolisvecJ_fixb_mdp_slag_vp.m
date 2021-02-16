% Calculate Coriolis joint torque vector for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:31
% EndTime: 2021-01-15 14:56:33
% DurationCPUTime: 0.51s
% Computational Cost: add. (304->112), mult. (733->159), div. (0->0), fcn. (374->4), ass. (0->57)
t100 = sin(qJ(3));
t103 = qJD(4) ^ 2;
t104 = qJD(3) ^ 2;
t142 = t100 * (t103 + t104);
t99 = sin(qJ(4));
t97 = t99 ^ 2;
t101 = cos(qJ(4));
t98 = t101 ^ 2;
t141 = (t97 - t98) * MDP(7);
t119 = MDP(12) + MDP(14);
t120 = MDP(11) + MDP(13);
t140 = t119 * t101 + t120 * t99;
t113 = (t97 + t98) * MDP(15);
t139 = pkin(4) * t97;
t138 = t101 * pkin(4);
t137 = qJ(5) + pkin(6);
t131 = qJD(4) * pkin(4);
t126 = t100 * qJD(2);
t89 = qJD(3) * pkin(6) + t126;
t112 = qJ(5) * qJD(3) + t89;
t125 = t101 * qJD(1);
t78 = -t112 * t99 - t125;
t75 = t78 + t131;
t136 = t75 - t78;
t122 = qJD(1) * qJD(4);
t128 = t99 * qJD(4);
t135 = t101 * t122 + t89 * t128;
t121 = qJD(3) * qJD(4);
t83 = t99 * pkin(4) * t121 + qJD(3) * t126;
t132 = qJD(3) * pkin(3);
t102 = cos(qJ(3));
t123 = t102 * qJD(2);
t95 = -pkin(3) - t138;
t129 = qJD(3) * t95;
t82 = qJD(5) - t123 + t129;
t130 = qJD(3) * t82;
t124 = t101 * qJD(4);
t118 = t101 * t104 * t99;
t117 = qJ(5) * t128;
t116 = qJ(5) * t124;
t114 = qJD(4) * t137;
t109 = -0.2e1 * t102 * t121;
t108 = qJD(5) + t123;
t96 = t99 * qJD(1);
t79 = t112 * t101 - t96;
t107 = t101 * t79 - t75 * t99;
t106 = -t108 - t82;
t90 = -t123 - t132;
t105 = qJD(3) * (-t90 - t123);
t93 = t99 * t122;
t86 = t137 * t101;
t85 = t137 * t99;
t81 = -t99 * qJD(5) - t101 * t114;
t80 = t101 * qJD(5) - t99 * t114;
t74 = -t89 * t124 + t93 + (-t108 * t99 - t116) * qJD(3);
t73 = (t108 * t101 - t117) * qJD(3) - t135;
t1 = [(-t107 * qJD(4) - t74 * t101 - t73 * t99) * MDP(16) + t140 * t103; t120 * (-t101 * t142 + t99 * t109) + t119 * (t101 * t109 + t99 * t142) + ((t107 * qJD(3) - t83) * MDP(16) + (-MDP(5) + t113) * t104) * t102 + (-t104 * MDP(4) + (t73 * t101 - t75 * t124 - t79 * t128 - t74 * t99 + t130) * MDP(16)) * t100; (t73 * t86 - t74 * t85 + t75 * t81 + t79 * t80 + t83 * t95) * MDP(16) + (t83 * MDP(14) + (-qJD(3) * t81 - t74) * MDP(15) + (pkin(6) * MDP(12) - MDP(9)) * t103) * t99 + (-t83 * MDP(13) + (qJD(3) * t80 + t73) * MDP(15) + (-pkin(6) * MDP(11) + MDP(8)) * t103) * t101 + ((-t100 * t82 - t102 * t107) * MDP(16) + (-t102 * t113 + (t101 * MDP(13) - t99 * MDP(14)) * t100) * qJD(3)) * qJD(2) + (t81 * MDP(13) - t80 * MDP(14) + (MDP(14) * t139 - 0.2e1 * t141) * qJD(3) + ((t90 - t132) * MDP(12) + (t82 + t129) * MDP(14) + (qJD(3) * t85 - t75) * MDP(15)) * t101 + (t90 * MDP(11) - t79 * MDP(15) + (pkin(4) * MDP(16) + MDP(13)) * t82 + (0.2e1 * t101 * MDP(6) - pkin(3) * MDP(11) + (t95 - t138) * MDP(13) - t86 * MDP(15)) * qJD(3)) * t99 + t140 * t123) * qJD(4); -MDP(6) * t118 + t104 * t141 + (-t96 * qJD(4) + t105 * t99 + t93) * MDP(11) + ((-t99 * t89 - t125) * qJD(4) + t101 * t105 + t135) * MDP(12) + (pkin(4) * t118 + t93 + (-t101 * t89 + t79) * qJD(4) + (t106 * t99 - t116) * qJD(3)) * MDP(13) + (-t104 * t139 + t78 * qJD(4) + (t101 * t106 + t117) * qJD(3) + t135) * MDP(14) + (-t131 + t136) * t101 * qJD(3) * MDP(15) + (t136 * t79 + (-t99 * t130 + t74) * pkin(4)) * MDP(16); t83 * MDP(16) - t104 * t113 + (-t107 * MDP(16) + 0.2e1 * (t99 * MDP(13) + t101 * MDP(14)) * qJD(4)) * qJD(3);];
tauc = t1;
