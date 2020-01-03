% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR8
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
%   see S4RPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:17
% EndTime: 2019-12-31 16:55:18
% DurationCPUTime: 0.40s
% Computational Cost: add. (304->90), mult. (695->138), div. (0->0), fcn. (406->4), ass. (0->53)
t113 = qJD(3) + qJD(4);
t117 = sin(qJ(3));
t119 = cos(qJ(3));
t155 = (t119 * MDP(12) - t117 * MDP(13)) * qJ(2) - t119 * t117 * MDP(7) + (t117 ^ 2 - t119 ^ 2) * MDP(8);
t154 = t117 * MDP(12) + t119 * MDP(13);
t152 = qJD(4) - t113;
t120 = -pkin(1) - pkin(5);
t150 = pkin(6) - t120;
t118 = cos(qJ(4));
t106 = t120 * qJD(1) + qJD(2);
t143 = qJD(1) * t117;
t90 = -pkin(6) * t143 + t117 * t106;
t149 = t118 * t90;
t147 = t118 * t119;
t116 = sin(qJ(4));
t132 = t116 * t143;
t146 = t113 * t132;
t142 = qJD(1) * t119;
t141 = qJD(3) * t117;
t140 = qJD(3) * t119;
t139 = qJD(4) * t116;
t134 = pkin(3) * t142;
t131 = t118 * t142;
t91 = -pkin(6) * t142 + t119 * t106;
t87 = qJD(3) * pkin(3) + t91;
t130 = -pkin(3) * t113 - t87;
t102 = t150 * t119;
t111 = t117 * pkin(3) + qJ(2);
t129 = pkin(6) * qJD(1) - t106;
t128 = -qJ(2) * MDP(6) - MDP(5);
t107 = pkin(3) * t140 + qJD(2);
t127 = t113 * t147;
t97 = t116 * t119 + t118 * t117;
t103 = t111 * qJD(1);
t88 = t129 * t141;
t89 = t129 * t140;
t92 = -t131 + t132;
t125 = t103 * t92 + t116 * t89 + t118 * t88;
t83 = t113 * t97;
t81 = t83 * qJD(1);
t93 = t97 * qJD(1);
t124 = -t92 * t93 * MDP(14) + (t93 * t113 - t81) * MDP(16) + (t146 + (-t131 - t92) * t113) * MDP(17) + (t92 ^ 2 - t93 ^ 2) * MDP(15);
t123 = t90 * t139 + (-t90 * t113 - t88) * t116 + t103 * t93;
t122 = qJD(1) ^ 2;
t121 = qJD(3) ^ 2;
t101 = t150 * t117;
t100 = t107 * qJD(1);
t98 = -t116 * t117 + t147;
t96 = qJD(3) * t102;
t95 = t150 * t141;
t84 = -t116 * t141 - t117 * t139 + t127;
t82 = qJD(1) * t127 - t146;
t1 = [(-t81 * t98 + t92 * t83) * MDP(14) + (t81 * t97 - t98 * t82 + t83 * t93 + t92 * t84) * MDP(15) + (t100 * t97 + t103 * t84 + t107 * t93 + t111 * t82) * MDP(19) + (t100 * t98 - t103 * t83 - t107 * t92 - t111 * t81) * MDP(20) + (-t83 * MDP(16) - t84 * MDP(17) + (t116 * t96 + t118 * t95) * MDP(19) + (-t116 * t95 + t118 * t96) * MDP(20) + ((t101 * t118 + t102 * t116) * MDP(19) + (-t101 * t116 + t102 * t118) * MDP(20)) * qJD(4)) * t113 + ((-MDP(13) * t120 - MDP(10)) * t119 + (-MDP(12) * t120 - MDP(9)) * t117) * t121 + (0.2e1 * (-t128 + t154) * qJD(2) + 0.2e1 * t155 * qJD(3)) * qJD(1); (-qJD(1) * t93 - t83 * t113) * MDP(19) + (qJD(1) * t92 - t84 * t113) * MDP(20) + t128 * t122 + t154 * (-t121 - t122); (-t93 * t134 - (-t116 * t91 - t149) * t113 + (t130 * t116 - t149) * qJD(4) + t125) * MDP(19) + (t92 * t134 + (t130 * qJD(4) + t91 * t113 + t89) * t118 + t123) * MDP(20) + t124 - t155 * t122; (t125 + t152 * (-t116 * t87 - t149)) * MDP(19) + ((-t152 * t87 + t89) * t118 + t123) * MDP(20) + t124;];
tauc = t1;
