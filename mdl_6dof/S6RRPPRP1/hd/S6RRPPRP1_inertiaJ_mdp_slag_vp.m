% Calculate joint inertia matrix for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:49
% EndTime: 2019-03-09 08:27:51
% DurationCPUTime: 0.71s
% Computational Cost: add. (1341->173), mult. (2444->250), div. (0->0), fcn. (2739->8), ass. (0->72)
t175 = -qJ(3) - pkin(7);
t130 = sin(pkin(10));
t132 = cos(pkin(10));
t134 = sin(qJ(5));
t168 = cos(qJ(5));
t174 = -t130 * t134 + t168 * t132;
t157 = t130 ^ 2 + t132 ^ 2;
t173 = t157 * MDP(15);
t117 = t168 * t130 + t134 * t132;
t150 = MDP(23) - MDP(26);
t151 = MDP(22) + MDP(24);
t153 = t132 * MDP(13);
t154 = t130 * MDP(14);
t172 = -t150 * t117 + t151 * t174 + t153 - t154;
t135 = sin(qJ(2));
t120 = t175 * t135;
t169 = cos(qJ(2));
t121 = t175 * t169;
t131 = sin(pkin(9));
t133 = cos(pkin(9));
t106 = -t133 * t120 - t121 * t131;
t171 = t106 ^ 2;
t170 = 2 * MDP(25);
t114 = t131 * t135 - t133 * t169;
t167 = pkin(5) * t114;
t123 = pkin(2) * t131 + qJ(4);
t166 = pkin(8) + t123;
t116 = t131 * t169 + t133 * t135;
t162 = t116 * t132;
t127 = -t169 * pkin(2) - pkin(1);
t104 = t114 * pkin(3) - t116 * qJ(4) + t127;
t108 = t120 * t131 - t121 * t133;
t93 = t132 * t104 - t108 * t130;
t90 = pkin(4) * t114 - pkin(8) * t162 + t93;
t163 = t116 * t130;
t94 = t130 * t104 + t132 * t108;
t92 = -pkin(8) * t163 + t94;
t86 = t134 * t90 + t168 * t92;
t99 = t174 * t116;
t165 = t174 * t99;
t164 = qJ(6) * t114;
t98 = t117 * t116;
t160 = t98 * MDP(20);
t159 = t99 * MDP(17);
t158 = t99 * MDP(19);
t156 = t114 * MDP(21);
t126 = -pkin(2) * t133 - pkin(3);
t155 = t126 * MDP(16);
t152 = 0.2e1 * t169;
t148 = t166 * t130;
t147 = t134 * t92 - t168 * t90;
t146 = -MDP(27) * pkin(5) - MDP(24);
t145 = t157 * MDP(16);
t96 = pkin(4) * t163 + t106;
t144 = -MDP(22) + t146;
t143 = t130 * t94 + t132 * t93;
t142 = -t130 * t93 + t132 * t94;
t141 = MDP(27) * qJ(6) - t150;
t140 = -t147 * MDP(22) - t86 * MDP(23);
t119 = -pkin(4) * t132 + t126;
t138 = MDP(13) * t130 + MDP(14) * t132;
t137 = t117 * MDP(19) + MDP(20) * t174;
t113 = t117 ^ 2;
t112 = t166 * t132;
t103 = t168 * t112 - t134 * t148;
t102 = t112 * t134 + t168 * t148;
t97 = -pkin(5) * t174 - qJ(6) * t117 + t119;
t95 = t117 * t98;
t87 = pkin(5) * t98 - qJ(6) * t99 + t96;
t84 = t147 - t167;
t83 = t164 + t86;
t1 = [MDP(1) + pkin(1) * MDP(9) * t152 + (t108 ^ 2 + t127 ^ 2 + t171) * MDP(12) + (t93 ^ 2 + t94 ^ 2 + t171) * MDP(16) + (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) * MDP(27) + (-0.2e1 * t98 * MDP(18) + t159) * t99 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t135 + MDP(5) * t152) * t135 + (t156 + 0.2e1 * t158 - 0.2e1 * t160) * t114 + (-t83 * t98 + t84 * t99) * t170 + 0.2e1 * (t98 * MDP(22) + t99 * MDP(23)) * t96 + 0.2e1 * (t98 * MDP(24) - t99 * MDP(26)) * t87 + 0.2e1 * (-t108 * MDP(11) + t93 * MDP(13) - t94 * MDP(14) - t84 * MDP(24) + t83 * MDP(26) + t140) * t114 + 0.2e1 * (-t143 * MDP(15) + (MDP(11) + t138) * t106) * t116; t135 * MDP(6) + t169 * MDP(7) + (-t106 * t132 + t126 * t163) * MDP(13) + (t106 * t130 + t126 * t162) * MDP(14) + t142 * MDP(15) + (t106 * t126 + t142 * t123) * MDP(16) + t117 * t159 + (-t95 + t165) * MDP(18) + (t119 * t98 - t174 * t96) * MDP(22) + (t117 * t96 + t119 * t99) * MDP(23) + (-t174 * t87 + t97 * t98) * MDP(24) + (t102 * t99 - t103 * t98 + t84 * t117 + t174 * t83) * MDP(25) + (-t117 * t87 - t97 * t99) * MDP(26) + (t102 * t84 + t103 * t83 + t87 * t97) * MDP(27) + (-t151 * t102 - t150 * t103 - t138 * t123 + t137) * t114 + (-t169 * MDP(10) - t135 * MDP(9)) * pkin(7) + ((-t114 * t131 - t116 * t133) * MDP(11) + (-t106 * t133 + t108 * t131) * MDP(12)) * pkin(2); MDP(8) + t113 * MDP(17) + (t102 ^ 2 + t103 ^ 2 + t97 ^ 2) * MDP(27) + (t131 ^ 2 + t133 ^ 2) * MDP(12) * pkin(2) ^ 2 + (-0.2e1 * t153 + 0.2e1 * t154 + t155) * t126 + (t102 * t117 + t103 * t174) * t170 + (t145 * t123 + 0.2e1 * t173) * t123 + 0.2e1 * (t119 * MDP(23) - t97 * MDP(26)) * t117 - 0.2e1 * (-MDP(18) * t117 + t119 * MDP(22) + t97 * MDP(24)) * t174; t127 * MDP(12) + t143 * MDP(16) + (-t95 - t165) * MDP(25) + (t117 * t83 - t174 * t84) * MDP(27) - t116 * t173 + t172 * t114; (-t102 * t174 + t103 * t117) * MDP(27); MDP(12) + t145 + (t174 ^ 2 + t113) * MDP(27); t106 * MDP(16) + t87 * MDP(27) + t138 * t116 + t150 * t99 + t151 * t98; t97 * MDP(27) + t155 - t172; 0; MDP(16) + MDP(27); t158 - t160 + t156 + (-t147 + 0.2e1 * t167) * MDP(24) + (-pkin(5) * t99 - qJ(6) * t98) * MDP(25) + (0.2e1 * t164 + t86) * MDP(26) + (-pkin(5) * t84 + qJ(6) * t83) * MDP(27) + t140; (-pkin(5) * t117 + qJ(6) * t174) * MDP(25) + t141 * t103 + t144 * t102 + t137; t141 * t117 - t144 * t174; 0; MDP(21) + 0.2e1 * pkin(5) * MDP(24) + 0.2e1 * qJ(6) * MDP(26) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(27); -t114 * MDP(24) + t99 * MDP(25) + t84 * MDP(27); t117 * MDP(25) + t102 * MDP(27); -t174 * MDP(27); 0; t146; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
