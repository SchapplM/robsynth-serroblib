% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:19
% EndTime: 2019-12-31 16:49:20
% DurationCPUTime: 0.40s
% Computational Cost: add. (316->82), mult. (781->137), div. (0->0), fcn. (480->6), ass. (0->53)
t108 = sin(pkin(7)) * pkin(1) + pkin(5);
t145 = pkin(6) + t108;
t111 = qJD(3) + qJD(4);
t146 = qJD(4) - t111;
t117 = sin(qJ(3));
t149 = MDP(5) * t117;
t116 = sin(qJ(4));
t118 = cos(qJ(4));
t119 = cos(qJ(3));
t99 = t116 * t119 + t118 * t117;
t148 = qJD(1) * t99;
t147 = (t117 ^ 2 - t119 ^ 2) * MDP(6);
t128 = t145 * qJD(1);
t88 = t119 * qJD(2) - t128 * t117;
t89 = t117 * qJD(2) + t128 * t119;
t144 = qJD(3) * pkin(3);
t143 = t118 * t89;
t120 = qJD(3) ^ 2;
t142 = t117 * t120;
t141 = t119 * t120;
t109 = -cos(pkin(7)) * pkin(1) - pkin(2);
t102 = qJD(1) * t109;
t139 = t117 * qJD(1);
t137 = t119 * qJD(1);
t136 = qJD(1) * qJD(3);
t135 = t117 * t144;
t134 = pkin(3) * t139;
t133 = t116 * t139;
t132 = t118 * t137;
t87 = t88 + t144;
t131 = -pkin(3) * t111 - t87;
t130 = t119 * t136;
t129 = qJD(3) * t145;
t98 = t116 * t117 - t118 * t119;
t100 = -t119 * pkin(3) + t109;
t125 = 0.2e1 * qJD(3) * t102;
t85 = t88 * qJD(3);
t86 = t89 * qJD(3);
t94 = -t116 * t137 - t118 * t139;
t95 = t100 * qJD(1);
t124 = -t116 * t85 - t118 * t86 + t95 * t94;
t79 = qJD(4) * t132 - t111 * t133 + t118 * t130;
t92 = -t132 + t133;
t123 = -t94 * t92 * MDP(12) + t79 * MDP(14) + (-t92 ^ 2 + t94 ^ 2) * MDP(13) + (t92 * MDP(14) + (-t148 - t94) * MDP(15)) * t111;
t122 = t95 * t92 + (t146 * t89 + t86) * t116;
t83 = t111 * t99;
t97 = t145 * t119;
t96 = t145 * t117;
t91 = t119 * t129;
t90 = t117 * t129;
t82 = t111 * t98;
t80 = t83 * qJD(1);
t1 = [0.2e1 * t130 * t149 - 0.2e1 * t136 * t147 + MDP(7) * t141 - MDP(8) * t142 + (-t108 * t141 + t117 * t125) * MDP(10) + (t108 * t142 + t119 * t125) * MDP(11) + (t79 * t99 + t94 * t82) * MDP(12) + (-t79 * t98 - t99 * t80 + t82 * t92 + t94 * t83) * MDP(13) + (t100 * t80 + t95 * t83 + (qJD(1) * t98 + t92) * t135) * MDP(17) + (t100 * t79 - t95 * t82 + (-t94 + t148) * t135) * MDP(18) + (-t82 * MDP(14) - t83 * MDP(15) + (t116 * t90 - t118 * t91 + (t116 * t96 - t118 * t97) * qJD(4)) * MDP(17) + (t116 * t91 + t118 * t90 - (-t116 * t97 - t118 * t96) * qJD(4)) * MDP(18)) * t111; (-t117 * MDP(10) - t119 * MDP(11)) * t120 + (-t83 * MDP(17) + t82 * MDP(18)) * t111; (-t92 * t134 - (-t116 * t88 - t143) * t111 + (t131 * t116 - t143) * qJD(4) + t124) * MDP(17) + (t94 * t134 + (t131 * qJD(4) + t88 * t111 - t85) * t118 + t122) * MDP(18) + t123 + (-t119 * t149 + t147) * qJD(1) ^ 2 + (-MDP(10) * t139 - MDP(11) * t137) * t102; (t124 + t146 * (-t116 * t87 - t143)) * MDP(17) + ((-t146 * t87 - t85) * t118 + t122) * MDP(18) + t123;];
tauc = t1;
