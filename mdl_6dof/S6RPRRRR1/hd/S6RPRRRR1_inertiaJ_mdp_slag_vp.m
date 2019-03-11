% Calculate joint inertia matrix for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR1_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:28
% EndTime: 2019-03-09 06:55:29
% DurationCPUTime: 0.43s
% Computational Cost: add. (798->129), mult. (1421->172), div. (0->0), fcn. (1594->10), ass. (0->77)
t130 = sin(qJ(6));
t134 = cos(qJ(6));
t165 = t130 * MDP(28) + t134 * MDP(29);
t163 = MDP(31) * t134;
t176 = -MDP(32) * t130 + t163;
t141 = 0.2e1 * t176;
t132 = sin(qJ(4));
t133 = sin(qJ(3));
t136 = cos(qJ(4));
t137 = cos(qJ(3));
t149 = t132 * t133 - t136 * t137;
t175 = t149 * MDP(17);
t129 = cos(pkin(11));
t158 = -t129 * pkin(1) - pkin(2);
t113 = -t137 * pkin(3) + t158;
t98 = t149 * pkin(4) + t113;
t174 = 0.2e1 * t98;
t173 = 0.2e1 * t113;
t172 = pkin(3) * t132;
t171 = pkin(5) * t130;
t135 = cos(qJ(5));
t170 = t135 * pkin(4);
t128 = sin(pkin(11));
t117 = pkin(1) * t128 + pkin(7);
t169 = pkin(8) + t117;
t131 = sin(qJ(5));
t112 = t132 * t137 + t133 * t136;
t110 = t169 * t133;
t111 = t169 * t137;
t155 = -t136 * t110 - t111 * t132;
t83 = -pkin(9) * t112 + t155;
t150 = t132 * t110 - t136 * t111;
t84 = -t149 * pkin(9) - t150;
t79 = t131 * t84 - t135 * t83;
t76 = t79 * t130;
t168 = t79 * t134;
t96 = t112 * t131 + t135 * t149;
t167 = MDP(30) * t96;
t166 = t130 * t134;
t97 = t135 * t112 - t131 * t149;
t94 = t97 * MDP(25);
t164 = MDP(25) * t131;
t121 = pkin(3) * t136 + pkin(4);
t114 = t135 * t121;
t103 = -t131 * t172 + t114;
t161 = t103 * MDP(24);
t104 = -t121 * t131 - t135 * t172;
t160 = t104 * MDP(25);
t159 = t136 * MDP(17);
t157 = MDP(27) * t166;
t126 = t130 ^ 2;
t156 = t126 * MDP(26) + MDP(23) + 0.2e1 * t157;
t154 = MDP(16) + t156;
t153 = -pkin(5) * t97 - pkin(10) * t96;
t101 = -pkin(5) - t103;
t102 = pkin(10) - t104;
t152 = t101 * t97 - t102 * t96;
t119 = pkin(4) * t131 + pkin(10);
t120 = -pkin(5) - t170;
t151 = -t119 * t96 + t120 * t97;
t148 = -t137 * MDP(10) + t133 * MDP(11);
t147 = MDP(28) * t134 - MDP(29) * t130;
t145 = -MDP(31) * t130 - MDP(32) * t134;
t144 = -t94 + (-MDP(24) - t176) * t96;
t143 = (MDP(24) * t135 - t164) * pkin(4);
t127 = t134 ^ 2;
t80 = t131 * t83 + t135 * t84;
t142 = -t79 * MDP(24) - t80 * MDP(25) + ((-t126 + t127) * MDP(27) + MDP(26) * t166 + MDP(21)) * t97 + (-MDP(22) + t165) * t96;
t140 = -t112 * MDP(18) + t144 - t175;
t139 = t112 * MDP(14) - t149 * MDP(15) + t155 * MDP(17) + t150 * MDP(18) + t142;
t125 = pkin(5) * t134;
t115 = t120 * t130;
t99 = t101 * t130;
t81 = t96 * pkin(5) - t97 * pkin(10) + t98;
t75 = t130 * t81 + t134 * t80;
t74 = -t130 * t80 + t134 * t81;
t1 = [MDP(1) + t173 * t175 + t94 * t174 + (t128 ^ 2 + t129 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(24) * t174 + t167) * t96 + 0.2e1 * (t74 * t96 + t97 * t76) * MDP(31) + 0.2e1 * (t97 * t168 - t75 * t96) * MDP(32) + (t127 * MDP(26) + MDP(19) - 0.2e1 * t157) * t97 ^ 2 + (MDP(12) * t112 - 0.2e1 * t149 * MDP(13) + MDP(18) * t173) * t112 + 0.2e1 * t148 * t158 + 0.2e1 * (-MDP(20) + t147) * t97 * t96 + (MDP(5) * t133 + 0.2e1 * t137 * MDP(6)) * t133; 0; MDP(4); (t152 * t130 - t168) * MDP(31) + (t152 * t134 + t76) * MDP(32) + (-MDP(10) * t133 - MDP(11) * t137) * t117 + t137 * MDP(8) + t133 * MDP(7) + t139; t140 - t148; MDP(9) - t101 * t141 + 0.2e1 * (-MDP(18) * t132 + t159) * pkin(3) + 0.2e1 * t161 + 0.2e1 * t160 + t154; (t151 * t130 - t168) * MDP(31) + (t151 * t134 + t76) * MDP(32) + t139; t140; (t114 + t170) * MDP(24) + (t115 + t99) * MDP(32) + (-t101 - t120) * t163 + (-pkin(4) - t121) * t164 + (t159 + (-MDP(24) * t131 - MDP(25) * t135 - MDP(18)) * t132) * pkin(3) + t154; -t120 * t141 + 0.2e1 * t143 + t154; (t153 * t130 - t168) * MDP(31) + (t153 * t134 + t76) * MDP(32) + t142; t144; t161 + t160 + (-t101 * t134 + t125) * MDP(31) + (t99 - t171) * MDP(32) + t156; (-t120 * t134 + t125) * MDP(31) + (t115 - t171) * MDP(32) + t143 + t156; pkin(5) * t141 + t156; t74 * MDP(31) - t75 * MDP(32) + t147 * t97 + t167; t145 * t97; t145 * t102 + t165; t145 * t119 + t165; t145 * pkin(10) + t165; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
