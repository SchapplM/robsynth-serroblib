% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:42:00
% EndTime: 2019-12-05 16:42:04
% DurationCPUTime: 0.68s
% Computational Cost: add. (538->121), mult. (979->172), div. (0->0), fcn. (434->4), ass. (0->64)
t136 = cos(qJ(3));
t175 = pkin(2) * qJD(3);
t151 = qJD(2) * t175;
t186 = qJD(1) * qJD(4) + t136 * t151;
t133 = sin(qJ(4));
t131 = t133 ^ 2;
t135 = cos(qJ(4));
t132 = t135 ^ 2;
t172 = t133 * t135;
t185 = MDP(8) * t172 - (t131 - t132) * MDP(9);
t130 = qJD(2) + qJD(3);
t134 = sin(qJ(3));
t176 = pkin(2) * qJD(2);
t153 = t134 * t176;
t116 = pkin(7) * t130 + t153;
t174 = t116 * t133;
t106 = qJD(1) * t135 - t174;
t184 = qJD(5) - t106;
t102 = -qJD(4) * pkin(4) + t184;
t107 = qJD(1) * t133 + t116 * t135;
t103 = qJD(4) * qJ(5) + t107;
t140 = t102 * t133 + t103 * t135;
t183 = (-MDP(7) + (t131 + t132) * MDP(16)) * t130 + t140 * MDP(18);
t163 = qJD(4) * t135;
t98 = t116 * t163 + t186 * t133;
t182 = -qJD(4) * t103 + t98;
t180 = t135 * MDP(13) + MDP(6);
t156 = MDP(13) + MDP(15);
t155 = MDP(14) - MDP(17);
t177 = pkin(2) * t136;
t118 = -pkin(4) * t135 - qJ(5) * t133 - pkin(3);
t173 = t118 * t130;
t152 = t136 * t176;
t164 = qJD(4) * t133;
t171 = t135 * t130 * t153 + t152 * t164;
t170 = t186 * t135;
t126 = pkin(2) * t134 + pkin(7);
t167 = MDP(18) * t126;
t162 = qJD(5) * t133;
t161 = t130 * MDP(14);
t137 = qJD(4) ^ 2;
t159 = t137 * MDP(11);
t108 = pkin(4) * t164 - qJ(5) * t163 - t162;
t104 = t134 * t175 + t108;
t141 = pkin(4) * t133 - qJ(5) * t135;
t143 = t134 * t151;
t97 = t143 + (qJD(4) * t141 - t162) * t130;
t149 = -t104 * t130 - t97;
t148 = -t108 * t130 - t97;
t146 = t106 + t174;
t145 = t156 * t137;
t144 = t155 * t137;
t139 = -t133 * t156 - t135 * t155;
t117 = -pkin(3) * t130 - t152;
t95 = (qJD(5) - t174) * qJD(4) + t170;
t138 = t137 * t135 * MDP(10) + (t117 * t163 + t133 * t143) * MDP(14) + (t102 * t163 + t98 * t133 + t95 * t135) * MDP(16) + 0.2e1 * t185 * qJD(4) * t130;
t129 = t130 ^ 2;
t127 = -pkin(3) - t177;
t113 = t118 - t177;
t110 = t117 * t164;
t109 = t141 * t130;
t101 = -t152 + t173;
t99 = t101 * t164;
t1 = [(qJD(4) * t140 + t95 * t133 - t98 * t135) * MDP(18) + t139 * t137; t110 * MDP(13) + t99 * MDP(15) + (t101 * t104 + t113 * t97) * MDP(18) + (t149 * MDP(15) + (t95 * MDP(18) - t145) * t126 + (t127 * t161 + (-t113 * t130 - t101) * MDP(17) + t102 * t167) * qJD(4)) * t135 + (-t159 + t149 * MDP(17) + (t98 * MDP(18) + t144) * t126 + ((MDP(13) * t127 + MDP(15) * t113) * t130 + (-MDP(16) - t167) * t103) * qJD(4)) * t133 + ((t133 * t161 + t180 * (-qJD(2) - t130)) * t134 + (-qJD(2) * MDP(7) + t139 * qJD(4) + t183) * t136) * t175 + t138; (t110 + t171) * MDP(13) + (t99 + t171) * MDP(15) + (t101 * t108 + t118 * t97) * MDP(18) + (t148 * MDP(15) + (-pkin(3) * t161 + (-t101 - t173) * MDP(17)) * qJD(4) + ((qJD(4) * t102 + t95) * MDP(18) - t145) * pkin(7)) * t135 + (-t159 + t148 * MDP(17) + (-t103 * MDP(16) + (-MDP(13) * pkin(3) + MDP(15) * t118) * t130) * qJD(4) + (t182 * MDP(18) + t144) * pkin(7)) * t133 + ((-t101 * MDP(18) - t180 * qJD(3) + (-t133 * t155 + MDP(6)) * t130) * t134 + (-qJD(3) * MDP(7) + t155 * t163 - t183) * t136) * t176 + t138; (qJ(5) * t95 - t101 * t109 - t102 * t107 + t184 * t103) * MDP(18) - t185 * t129 + (t146 * MDP(14) + (0.2e1 * qJD(5) - t146) * MDP(17) + t156 * t107) * qJD(4) + ((-t117 * MDP(13) - t101 * MDP(15) + t109 * MDP(17)) * t133 + (-t117 * MDP(14) + t109 * MDP(15) + t101 * MDP(17)) * t135) * t130 + (-MDP(18) * pkin(4) - t156) * t98 - t155 * t170; -t129 * MDP(15) * t172 + (-t129 * t131 - t137) * MDP(17) + (t101 * t130 * t133 + t182) * MDP(18);];
tauc = t1;
