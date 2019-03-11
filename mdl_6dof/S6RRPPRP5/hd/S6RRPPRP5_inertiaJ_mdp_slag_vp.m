% Calculate joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:59
% EndTime: 2019-03-09 08:44:01
% DurationCPUTime: 0.68s
% Computational Cost: add. (848->172), mult. (1418->233), div. (0->0), fcn. (1367->6), ass. (0->64)
t174 = pkin(3) + pkin(7);
t133 = sin(pkin(9));
t134 = cos(pkin(9));
t136 = sin(qJ(5));
t168 = cos(qJ(5));
t111 = t136 * t133 - t168 * t134;
t138 = cos(qJ(2));
t103 = t111 * t138;
t143 = -t168 * t133 - t136 * t134;
t104 = t143 * t138;
t135 = -pkin(2) - qJ(4);
t137 = sin(qJ(2));
t155 = -t137 * qJ(3) - pkin(1);
t109 = t135 * t138 + t155;
t118 = t174 * t137;
t115 = t134 * t118;
t95 = t137 * pkin(4) + t115 + (pkin(8) * t138 - t109) * t133;
t100 = t134 * t109 + t133 * t118;
t165 = t134 * t138;
t97 = -pkin(8) * t165 + t100;
t154 = t136 * t97 - t168 * t95;
t90 = t136 * t95 + t168 * t97;
t173 = t104 * MDP(21) + t103 * MDP(22) - t154 * MDP(24) - t90 * MDP(25);
t170 = -2 * MDP(20);
t169 = 2 * MDP(28);
t167 = t137 * pkin(5);
t166 = -pkin(8) + t135;
t163 = t137 * qJ(6);
t99 = -t133 * t109 + t115;
t91 = t100 * t133 + t99 * t134;
t162 = t91 * MDP(18);
t123 = t133 * pkin(4) + qJ(3);
t119 = t174 * t138;
t121 = t133 ^ 2 + t134 ^ 2;
t161 = t104 * MDP(19);
t160 = MDP(24) + MDP(26);
t159 = -MDP(25) + MDP(28);
t107 = pkin(4) * t165 + t119;
t157 = t111 ^ 2 + t143 ^ 2;
t156 = t166 * t134;
t153 = -pkin(2) * MDP(14) + MDP(12);
t152 = -MDP(29) * pkin(5) - MDP(26);
t151 = -MDP(24) + t152;
t87 = t163 + t90;
t88 = t154 - t167;
t150 = t88 * t111 - t143 * t87;
t149 = MDP(29) * qJ(6) + t159;
t147 = t134 * MDP(15) - t133 * MDP(16);
t146 = t133 * MDP(15) + t134 * MDP(16);
t144 = -t111 * MDP(21) + MDP(22) * t143;
t142 = MDP(11) + t147;
t141 = qJ(3) * MDP(18) + t146;
t140 = pkin(7) ^ 2;
t139 = qJ(3) ^ 2;
t132 = t138 ^ 2;
t131 = t137 ^ 2;
t117 = -t138 * pkin(2) + t155;
t116 = t166 * t133;
t110 = t121 * t135;
t102 = t168 * t116 + t136 * t156;
t101 = t136 * t116 - t168 * t156;
t98 = -pkin(5) * t143 + t111 * qJ(6) + t123;
t92 = -t103 * pkin(5) - t104 * qJ(6) + t107;
t1 = [(t87 ^ 2 + t88 ^ 2 + t92 ^ 2) * MDP(29) + (t117 ^ 2 + t132 * t140) * MDP(14) + (t100 ^ 2 + t119 ^ 2 + t99 ^ 2) * MDP(18) + MDP(1) + (t140 * MDP(14) + MDP(23) + MDP(4)) * t131 - (t103 * t170 - t161) * t104 - 0.2e1 * (-t107 * MDP(25) - t88 * MDP(27) + t92 * MDP(28)) * t104 + 0.2e1 * (-t107 * MDP(24) - t92 * MDP(26) + t87 * MDP(27)) * t103 + 0.2e1 * (t131 + t132) * MDP(11) * pkin(7) + 0.2e1 * (-pkin(1) * MDP(10) - t117 * MDP(13) + t99 * MDP(15) - t100 * MDP(16) - t88 * MDP(26) + t87 * MDP(28) + t138 * MDP(5) + t173) * t137 + 0.2e1 * ((-t100 * t134 + t133 * t99) * MDP(17) + t147 * t119 + t117 * MDP(12) + pkin(1) * MDP(9)) * t138; -t91 * MDP(17) - t111 * t161 + (-t111 * t103 + t104 * t143) * MDP(20) + (-t123 * t103 - t107 * t143) * MDP(24) + (t123 * t104 - t107 * t111) * MDP(25) + (-t98 * t103 - t143 * t92) * MDP(26) + (t101 * t104 + t102 * t103 - t150) * MDP(27) + (-t98 * t104 + t92 * t111) * MDP(28) + (t88 * t101 + t87 * t102 + t92 * t98) * MDP(29) + t135 * t162 + t141 * t119 + (t142 * qJ(3) + MDP(7)) * t138 + (-pkin(2) * MDP(11) - t160 * t101 + t159 * t102 + t147 * t135 + MDP(6) + t144) * t137 + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t138 + (-MDP(9) + t153) * t137) * pkin(7); MDP(8) - 0.2e1 * pkin(2) * MDP(12) + (pkin(2) ^ 2 + t139) * MDP(14) - 0.2e1 * t110 * MDP(17) + (t121 * t135 ^ 2 + t139) * MDP(18) + (t101 ^ 2 + t102 ^ 2 + t98 ^ 2) * MDP(29) + (MDP(19) * t111 - 0.2e1 * t123 * MDP(25) - 0.2e1 * t101 * MDP(27) + t143 * t170 + t98 * t169) * t111 - 0.2e1 * (t123 * MDP(24) + t98 * MDP(26) - t102 * MDP(27)) * t143 + 0.2e1 * (MDP(13) + t146) * qJ(3); t162 + (-t103 * t143 + t111 * t104) * MDP(27) + t150 * MDP(29) + (pkin(7) * MDP(14) - t160 * t111 - t143 * t159 + t142) * t137; -t121 * MDP(17) + t110 * MDP(18) - t157 * MDP(27) + (t101 * t111 - t102 * t143) * MDP(29) + t153; t121 * MDP(18) + t157 * MDP(29) + MDP(14); t119 * MDP(18) + t92 * MDP(29) - t160 * t103 - t159 * t104 + t147 * t138; t98 * MDP(29) + t159 * t111 - t143 * t160 + t141; 0; MDP(18) + MDP(29); t137 * MDP(23) + (-t154 + 0.2e1 * t167) * MDP(26) + (-pkin(5) * t104 + t103 * qJ(6)) * MDP(27) + (0.2e1 * t163 + t90) * MDP(28) + (-pkin(5) * t88 + t87 * qJ(6)) * MDP(29) + t173; (t111 * pkin(5) + qJ(6) * t143) * MDP(27) + t149 * t102 + t151 * t101 + t144; t151 * t111 - t143 * t149; 0; MDP(23) + 0.2e1 * pkin(5) * MDP(26) + qJ(6) * t169 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29); -t137 * MDP(26) + t104 * MDP(27) + t88 * MDP(29); -t111 * MDP(27) + t101 * MDP(29); t111 * MDP(29); 0; t152; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
