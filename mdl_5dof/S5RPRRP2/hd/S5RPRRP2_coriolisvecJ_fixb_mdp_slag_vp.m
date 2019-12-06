% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:46
% EndTime: 2019-12-05 18:01:48
% DurationCPUTime: 0.52s
% Computational Cost: add. (657->117), mult. (1317->167), div. (0->0), fcn. (688->6), ass. (0->72)
t184 = qJ(5) + pkin(7);
t142 = sin(qJ(4));
t138 = t142 ^ 2;
t144 = cos(qJ(4));
t139 = t144 ^ 2;
t183 = (t138 - t139) * MDP(9);
t131 = cos(pkin(8)) * pkin(1) + pkin(2);
t143 = sin(qJ(3));
t145 = cos(qJ(3));
t180 = pkin(1) * sin(pkin(8));
t182 = t145 * t131 - t143 * t180;
t125 = t131 * qJD(1);
t162 = qJD(1) * t180;
t113 = t143 * t125 + t145 * t162;
t137 = qJD(1) + qJD(3);
t156 = t184 * t137 + t113;
t99 = t144 * qJD(2) - t156 * t142;
t100 = t142 * qJD(2) + t156 * t144;
t181 = qJD(4) * t100;
t179 = t137 * pkin(3);
t178 = t144 * pkin(4);
t175 = qJD(4) * pkin(4);
t98 = t99 + t175;
t176 = t98 - t99;
t174 = t113 * t137;
t146 = qJD(4) ^ 2;
t173 = t142 * t146;
t172 = t144 * t146;
t148 = t143 * t131 + t145 * t180;
t118 = pkin(7) + t148;
t170 = -qJ(5) - t118;
t112 = t145 * t125 - t143 * t162;
t106 = -t112 - t179;
t152 = qJD(3) * t162;
t166 = qJD(3) * t125;
t111 = t143 * t166 + t145 * t152;
t163 = t144 * qJD(4);
t169 = t106 * t163 + t111 * t142;
t168 = -t138 - t139;
t164 = t142 * qJD(4);
t161 = pkin(4) * t164;
t159 = t137 * t164;
t102 = pkin(4) * t159 + t111;
t160 = -pkin(3) - t178;
t110 = -t143 * t152 + t145 * t166;
t153 = qJD(5) * t137 + t110;
t94 = t99 * qJD(4) + t153 * t144;
t95 = -t153 * t142 - t181;
t158 = -t95 * t142 + t94 * t144;
t157 = qJD(4) * t184;
t155 = -t106 * t137 - t110;
t154 = qJD(4) * t170;
t117 = -pkin(3) - t182;
t151 = pkin(7) * t146 - t174;
t150 = qJD(4) * (t112 - t179);
t147 = 0.2e1 * t144 * MDP(8) * t159 - 0.2e1 * t137 * qJD(4) * t183 + MDP(10) * t172 - MDP(11) * t173;
t116 = t148 * qJD(3);
t136 = t137 ^ 2;
t135 = t144 * qJ(5);
t133 = t144 * qJD(5);
t127 = t144 * pkin(7) + t135;
t126 = t184 * t142;
t120 = -t142 * qJD(5) - t144 * t157;
t119 = -t142 * t157 + t133;
t115 = t182 * qJD(3);
t109 = t144 * t118 + t135;
t108 = t170 * t142;
t103 = t106 * t164;
t101 = t160 * t137 + qJD(5) - t112;
t97 = (-qJD(5) - t115) * t142 + t144 * t154;
t96 = t144 * t115 + t142 * t154 + t133;
t1 = [-t111 * MDP(6) - t110 * MDP(7) + (-t111 * t144 - t115 * t164 - t118 * t172 + t103) * MDP(13) + (-t115 * t163 + t118 * t173 + t169) * MDP(14) + (-t100 * t164 - t98 * t163 + t158) * MDP(15) + (t94 * t109 + t100 * t96 + t95 * t108 + t98 * t97 + t102 * (t117 - t178) + t101 * (t116 + t161)) * MDP(16) + (-t116 * MDP(6) - t115 * MDP(7) + (-t116 * t144 + t117 * t164) * MDP(13) + (t116 * t142 + t117 * t163) * MDP(14) + (-t108 * t163 - t109 * t164 - t142 * t97 + t144 * t96) * MDP(15)) * t137 + t147; (-t146 * MDP(14) + (t95 + t181) * MDP(16)) * t144 + (-t146 * MDP(13) + (-qJD(4) * t98 + t94) * MDP(16)) * t142; (-t111 + t174) * MDP(6) + (t112 * t137 - t110) * MDP(7) + (t103 + t142 * t150 + (-t111 - t151) * t144) * MDP(13) + (t151 * t142 + t144 * t150 + t169) * MDP(14) + ((-t100 * t142 - t144 * t98) * qJD(4) + (t119 * t144 - t120 * t142 + t168 * t112 + (t126 * t144 - t127 * t142) * qJD(4)) * t137 + t158) * MDP(15) + (t94 * t127 - t95 * t126 + t102 * t160 + (t142 * t112 + t120) * t98 + (-t113 + t161) * t101 + (-t144 * t112 + t119) * t100) * MDP(16) + t147; t136 * t183 + t155 * t142 * MDP(13) + (t176 * t100 + (-t101 * t137 * t142 + t95) * pkin(4)) * MDP(16) + (-t142 * t136 * MDP(8) + t155 * MDP(14) + (-t175 + t176) * t137 * MDP(15)) * t144; ((-t100 * t144 + t142 * t98) * t137 + t102) * MDP(16) + t168 * MDP(15) * t136;];
tauc = t1;
