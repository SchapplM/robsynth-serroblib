% Calculate joint inertia matrix for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:43
% EndTime: 2019-03-09 08:39:46
% DurationCPUTime: 0.80s
% Computational Cost: add. (804->193), mult. (1476->258), div. (0->0), fcn. (1398->6), ass. (0->71)
t175 = pkin(7) * MDP(14);
t125 = sin(pkin(9));
t126 = cos(pkin(9));
t148 = MDP(12) - MDP(17);
t149 = MDP(11) + MDP(15);
t174 = t149 * t125 + t148 * t126;
t127 = sin(qJ(5));
t129 = cos(qJ(5));
t173 = -t125 * t129 + t127 * t126;
t172 = MDP(14) + MDP(18);
t128 = sin(qJ(2));
t171 = 0.2e1 * t128;
t170 = -2 * MDP(20);
t169 = pkin(3) + pkin(4);
t130 = cos(qJ(2));
t168 = pkin(5) * t130;
t167 = pkin(7) * t130;
t166 = pkin(8) * t128;
t165 = -pkin(8) + qJ(3);
t111 = -pkin(2) * t130 - qJ(3) * t128 - pkin(1);
t116 = t125 * t167;
t121 = t130 * pkin(3);
t91 = pkin(4) * t130 + t116 + t121 + (-t111 - t166) * t126;
t100 = t125 * t111 + t126 * t167;
t97 = -qJ(4) * t130 + t100;
t93 = t125 * t166 + t97;
t85 = t127 * t91 + t129 * t93;
t164 = pkin(2) * MDP(14);
t163 = qJ(6) * t130;
t161 = t126 * t128;
t144 = qJ(4) * t125 + pkin(2);
t109 = -pkin(3) * t126 - t144;
t158 = MDP(18) * t109;
t102 = t173 * t128;
t157 = t102 * MDP(22);
t107 = t125 * t127 + t126 * t129;
t103 = t107 * t128;
t156 = t103 * MDP(19);
t155 = t103 * MDP(21);
t154 = t125 * MDP(12);
t153 = t125 * MDP(17);
t152 = t126 * MDP(11);
t151 = t126 * MDP(15);
t150 = t130 * MDP(23);
t147 = MDP(24) + MDP(26);
t146 = -MDP(25) + MDP(28);
t145 = t165 * t125;
t143 = t127 * t93 - t129 * t91;
t142 = -MDP(29) * pkin(5) - MDP(26);
t99 = t111 * t126 - t116;
t139 = MDP(24) - t142;
t98 = t121 - t99;
t138 = t125 * t98 + t126 * t97;
t137 = t100 * t126 - t125 * t99;
t136 = MDP(29) * qJ(6) + t146;
t135 = -t143 * MDP(24) - t85 * MDP(25);
t134 = t125 * MDP(11) + MDP(12) * t126;
t133 = t125 * MDP(15) - MDP(17) * t126;
t132 = -MDP(21) * t173 - t107 * MDP(22);
t104 = t169 * t126 + t144;
t115 = qJ(4) * t161;
t96 = t115 + (-t169 * t125 - pkin(7)) * t128;
t112 = t165 * t126;
t101 = -t115 + (pkin(3) * t125 + pkin(7)) * t128;
t95 = t129 * t112 + t127 * t145;
t94 = t112 * t127 - t129 * t145;
t88 = pkin(5) * t107 + qJ(6) * t173 + t104;
t86 = -pkin(5) * t102 + qJ(6) * t103 - t96;
t83 = t143 - t168;
t82 = t163 + t85;
t1 = [(t100 ^ 2 + t99 ^ 2) * MDP(14) + (t101 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(18) + (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) * MDP(29) + MDP(1) + (t102 * t170 + t156) * t103 + (MDP(5) * t171 + 0.2e1 * pkin(1) * MDP(9) + t150 + 0.2e1 * t155 - 0.2e1 * t157) * t130 + 0.2e1 * (t96 * MDP(25) + t83 * MDP(27) + t86 * MDP(28)) * t103 + 0.2e1 * (t96 * MDP(24) - t86 * MDP(26) - t82 * MDP(27)) * t102 + 0.2e1 * (-t99 * MDP(11) + t100 * MDP(12) + t98 * MDP(15) - t97 * MDP(17) - t83 * MDP(26) + t82 * MDP(28) + t135) * t130 + ((-t100 * t125 - t126 * t99) * MDP(13) + (-t125 * t97 + t126 * t98) * MDP(16) + t133 * t101) * t171 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(4) + (0.2e1 * t134 + t175) * pkin(7)) * t128) * t128; t137 * MDP(13) + t138 * MDP(16) - t173 * t156 + (t102 * t173 - t103 * t107) * MDP(20) + (t102 * t104 + t107 * t96) * MDP(24) + (t103 * t104 - t173 * t96) * MDP(25) + (t102 * t88 - t107 * t86) * MDP(26) + (-t102 * t95 + t103 * t94 - t107 * t82 - t173 * t83) * MDP(27) + (-t103 * t88 - t173 * t86) * MDP(28) + (t82 * t95 + t83 * t94 - t86 * t88) * MDP(29) + (-t151 - t153 + t158) * t101 + (-pkin(7) * MDP(10) + t146 * t95 - t147 * t94 + MDP(7) + t132) * t130 + (MDP(6) + t133 * t109 - t134 * pkin(2) + (-MDP(9) - t152 + t154 - t164) * pkin(7)) * t128 + (t137 * MDP(14) + t138 * MDP(18) + t130 * t174) * qJ(3); MDP(8) + (t88 ^ 2 + t94 ^ 2 + t95 ^ 2) * MDP(29) + (-0.2e1 * t151 - 0.2e1 * t153 + t158) * t109 + 0.2e1 * (t104 * MDP(24) + t88 * MDP(26)) * t107 + (0.2e1 * t152 - 0.2e1 * t154 + t164) * pkin(2) - (-MDP(19) * t173 + 0.2e1 * MDP(25) * t104 - 0.2e1 * MDP(28) * t88 + t107 * t170) * t173 + 0.2e1 * (-t107 * t95 - t173 * t94) * MDP(27) + (t172 * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * MDP(16)) * (t125 ^ 2 + t126 ^ 2) * qJ(3); t101 * MDP(18) + t86 * MDP(29) + t146 * t103 - t147 * t102 + (t174 + t175) * t128; -MDP(29) * t88 - t147 * t107 + t148 * t125 - t149 * t126 - t146 * t173 + t158 - t164; MDP(29) + t172; MDP(16) * t161 + t98 * MDP(18) + (-t102 * t127 - t103 * t129) * MDP(27) + (t127 * t82 - t129 * t83) * MDP(29) + (t146 * t127 + t147 * t129 + MDP(15)) * t130; (-t107 * t127 + t129 * t173) * MDP(27) + (t127 * t95 - t129 * t94) * MDP(29) + (MDP(18) * qJ(3) + MDP(16)) * t125; 0; MDP(18) + (t127 ^ 2 + t129 ^ 2) * MDP(29); t155 - t157 + t150 + (-t143 + 0.2e1 * t168) * MDP(26) + (-pkin(5) * t103 - qJ(6) * t102) * MDP(27) + (0.2e1 * t163 + t85) * MDP(28) + (-pkin(5) * t83 + qJ(6) * t82) * MDP(29) + t135; (pkin(5) * t173 - qJ(6) * t107) * MDP(27) + t136 * t95 - t139 * t94 + t132; 0; t136 * t127 + t139 * t129; MDP(23) + 0.2e1 * pkin(5) * MDP(26) + 0.2e1 * qJ(6) * MDP(28) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29); -t130 * MDP(26) + t103 * MDP(27) + t83 * MDP(29); -MDP(27) * t173 + t94 * MDP(29); 0; -t129 * MDP(29); t142; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
