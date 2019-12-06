% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:55
% EndTime: 2019-12-05 16:06:59
% DurationCPUTime: 0.73s
% Computational Cost: add. (772->156), mult. (1925->212), div. (0->0), fcn. (1195->4), ass. (0->77)
t199 = MDP(12) + MDP(15);
t178 = (qJD(2) * qJD(3));
t198 = -2 * t178;
t158 = sin(qJ(3));
t197 = MDP(5) * t158;
t159 = cos(qJ(3));
t196 = MDP(6) * (t158 ^ 2 - t159 ^ 2);
t194 = t158 * MDP(10) + t159 * MDP(11);
t156 = sin(pkin(8));
t157 = cos(pkin(8));
t142 = t156 * t159 + t157 * t158;
t138 = t142 * qJD(2);
t134 = t138 ^ 2;
t192 = -qJ(4) - pkin(6);
t188 = t157 * t159;
t141 = t156 * t158 - t188;
t170 = qJD(3) * t192;
t131 = qJD(4) * t159 + t158 * t170;
t177 = qJD(3) * qJD(1);
t120 = qJD(2) * t131 + t159 * t177;
t164 = -qJD(4) * t158 + t159 * t170;
t162 = qJD(2) * t164 - t158 * t177;
t99 = t120 * t156 - t157 * t162;
t191 = t141 * t99;
t145 = t192 * t159;
t173 = t192 * t158;
t121 = -t145 * t156 - t157 * t173;
t190 = t99 * t121;
t132 = qJD(1) * t158 - qJD(2) * t145;
t189 = t132 * t156;
t124 = t157 * t132;
t160 = qJD(3) ^ 2;
t187 = t158 * t160;
t186 = t159 * t160;
t100 = t157 * t120 + t156 * t162;
t130 = t159 * qJD(1) + qJD(2) * t173;
t127 = qJD(3) * pkin(3) + t130;
t108 = t156 * t127 + t124;
t184 = qJD(2) * t159;
t183 = qJD(3) * t158;
t181 = t158 * qJD(2);
t112 = t130 * t157 - t189;
t179 = qJD(5) - t112;
t176 = pkin(3) * t183;
t175 = -pkin(3) * t159 - pkin(2);
t174 = t157 * t184;
t172 = t158 * t178;
t171 = t159 * t178;
t169 = pkin(2) * t198;
t168 = t175 * qJD(2);
t107 = t127 * t157 - t189;
t137 = t142 * qJD(3);
t128 = qJD(2) * t137;
t146 = t156 * t172;
t129 = t157 * t171 - t146;
t149 = pkin(3) * t172;
t167 = pkin(4) * t128 - qJ(5) * t129 + t149;
t135 = t156 * t181 - t174;
t144 = qJD(4) + t168;
t106 = pkin(4) * t135 - qJ(5) * t138 + t144;
t166 = t106 * t138 + t99;
t111 = t131 * t156 - t157 * t164;
t113 = t157 * t131 + t156 * t164;
t122 = -t157 * t145 + t156 * t173;
t163 = t111 * t138 - t113 * t135 + t121 * t129 - t122 * t128 + t142 * t99;
t152 = -pkin(3) * t157 - pkin(4);
t150 = pkin(3) * t156 + qJ(5);
t140 = qJD(3) * t188 - t156 * t183;
t114 = pkin(4) * t141 - qJ(5) * t142 + t175;
t110 = t130 * t156 + t124;
t109 = pkin(3) * t181 + pkin(4) * t138 + qJ(5) * t135;
t104 = qJD(3) * qJ(5) + t108;
t103 = -qJD(3) * pkin(4) + qJD(5) - t107;
t102 = pkin(4) * t137 - qJ(5) * t140 - qJD(5) * t142 + t176;
t98 = qJD(3) * qJD(5) + t100;
t97 = -qJD(5) * t138 + t167;
t1 = [(t100 * t142 - t107 * t137 + t108 * t140 + t191) * MDP(13) + (t103 * t137 + t104 * t140 + t142 * t98 + t191) * MDP(17) - t194 * t160 + (-MDP(14) * t137 + MDP(16) * t140) * qJD(3) + t199 * (-t142 * t128 + t129 * t141 - t140 * t135 + t137 * t138); 0.2e1 * t171 * t197 + t196 * t198 + MDP(7) * t186 - MDP(8) * t187 + (-pkin(6) * t186 + t158 * t169) * MDP(10) + (pkin(6) * t187 + t159 * t169) * MDP(11) + (-t100 * t141 - t107 * t140 - t108 * t137 + t163) * MDP(12) + (t100 * t122 - t107 * t111 + t108 * t113 + t190 + (t144 + t168) * t176) * MDP(13) + (-qJD(3) * t111 + t102 * t135 + t106 * t137 + t114 * t128 + t141 * t97) * MDP(14) + (t103 * t140 - t104 * t137 - t141 * t98 + t163) * MDP(15) + (qJD(3) * t113 - t102 * t138 - t106 * t140 - t114 * t129 - t142 * t97) * MDP(16) + (t102 * t106 + t103 * t111 + t104 * t113 + t114 * t97 + t122 * t98 + t190) * MDP(17); ((t108 - t110) * t138 + (-t107 + t112) * t135 + (-t128 * t156 - t129 * t157) * pkin(3)) * MDP(12) + (t107 * t110 - t108 * t112 + (t100 * t156 - t144 * t181 - t157 * t99) * pkin(3)) * MDP(13) + (qJD(3) * t110 - t109 * t135 - t166) * MDP(14) + (-t128 * t150 + t152 * t129 + (t104 - t110) * t138 + (t103 - t179) * t135) * MDP(15) + (-t106 * t135 + t109 * t138 + (0.2e1 * qJD(5) - t112) * qJD(3) + t100) * MDP(16) + (-t103 * t110 + t104 * t179 - t106 * t109 + t150 * t98 + t152 * t99) * MDP(17) + (pkin(2) * t194 - t159 * t197 + t196) * qJD(2) ^ 2; (t107 * t138 + t108 * t135 + t149) * MDP(13) + t146 * MDP(16) + (t104 * t135 + (-qJD(5) - t103) * t138 + t167) * MDP(17) + ((t156 * t184 + t157 * t181 + t138) * MDP(14) + (t135 - t174) * MDP(16)) * qJD(3) + t199 * (-t135 ^ 2 - t134); t138 * t135 * MDP(14) + (-t146 + (t135 + t174) * qJD(3)) * MDP(15) + (-t134 - t160) * MDP(16) + (-qJD(3) * t104 + t166) * MDP(17);];
tauc = t1;
