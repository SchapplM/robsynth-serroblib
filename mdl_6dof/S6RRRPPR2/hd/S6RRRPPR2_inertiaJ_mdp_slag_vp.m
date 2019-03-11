% Calculate joint inertia matrix for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:44
% EndTime: 2019-03-09 15:26:46
% DurationCPUTime: 0.54s
% Computational Cost: add. (934->147), mult. (1682->201), div. (0->0), fcn. (1866->8), ass. (0->72)
t134 = sin(qJ(6));
t137 = cos(qJ(6));
t176 = MDP(29) * t134 + MDP(30) * t137;
t177 = MDP(22) + t176;
t135 = sin(qJ(3));
t136 = sin(qJ(2));
t138 = cos(qJ(3));
t139 = cos(qJ(2));
t114 = t135 * t136 - t138 * t139;
t115 = t135 * t139 + t136 * t138;
t132 = sin(pkin(10));
t133 = cos(pkin(10));
t100 = t133 * t114 + t115 * t132;
t101 = -t114 * t132 + t115 * t133;
t127 = -pkin(2) * t139 - pkin(1);
t104 = pkin(3) * t114 + t127;
t142 = -qJ(5) * t101 + t104;
t89 = pkin(4) * t100 + t142;
t175 = -0.2e1 * t89;
t174 = 0.2e1 * t127;
t173 = 0.2e1 * t139;
t172 = 2 * MDP(21);
t171 = pkin(7) + pkin(8);
t170 = pkin(2) * t135;
t117 = t171 * t136;
t118 = t171 * t139;
t153 = -t138 * t117 - t118 * t135;
t145 = -qJ(4) * t115 + t153;
t148 = t117 * t135 - t118 * t138;
t97 = -qJ(4) * t114 - t148;
t92 = t132 * t145 + t133 * t97;
t86 = -t100 * pkin(5) + t92;
t82 = t86 * t134;
t83 = t86 * t137;
t126 = pkin(2) * t138 + pkin(3);
t167 = t126 * t132;
t110 = t133 * t170 + t167;
t106 = qJ(5) + t110;
t169 = t100 * t106;
t122 = pkin(3) * t132 + qJ(5);
t168 = t100 * t122;
t166 = t134 * t137;
t109 = t126 * t133 - t132 * t170;
t108 = -pkin(4) - t109;
t164 = MDP(23) * t108;
t124 = -pkin(3) * t133 - pkin(4);
t163 = MDP(23) * t124;
t162 = MDP(28) * t101;
t159 = t100 * MDP(21);
t158 = t101 * MDP(27);
t129 = t137 * MDP(26);
t157 = t138 * MDP(16);
t90 = t132 * t97 - t133 * t145;
t156 = t90 ^ 2 + t92 ^ 2;
t155 = MDP(25) * t166;
t131 = t137 ^ 2;
t154 = t131 * MDP(24) + MDP(15) - 0.2e1 * t155;
t152 = -MDP(27) * t134 + t129;
t105 = -pkin(9) + t108;
t150 = -t101 * t105 + t169;
t121 = -pkin(9) + t124;
t149 = -t101 * t121 + t168;
t147 = MDP(29) * t137 - t134 * MDP(30);
t144 = 0.2e1 * t177;
t143 = (MDP(26) * t134 + MDP(27) * t137) * t100;
t130 = t134 ^ 2;
t141 = t115 * MDP(13) - t114 * MDP(14) + t153 * MDP(16) + t148 * MDP(17) + t90 * MDP(21) + t92 * MDP(22) + t83 * MDP(30) + t101 * t129 + ((-t130 + t131) * MDP(25) + MDP(24) * t166) * t100;
t85 = pkin(5) * t101 + t90;
t84 = (pkin(4) + pkin(9)) * t100 + t142;
t81 = t134 * t85 + t137 * t84;
t80 = -t134 * t84 + t137 * t85;
t1 = [MDP(1) + pkin(1) * MDP(9) * t173 + t114 * MDP(16) * t174 + (t104 ^ 2 + t156) * MDP(19) + t159 * t175 + (t89 ^ 2 + t156) * MDP(23) + (t130 * MDP(24) + 0.2e1 * t155) * t100 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t136 + MDP(5) * t173) * t136 + (MDP(11) * t115 - 0.2e1 * t114 * MDP(12) + MDP(17) * t174) * t115 + (MDP(22) * t175 + 0.2e1 * t143 + t162) * t101 + 0.2e1 * (-t100 * t83 + t101 * t80) * MDP(29) + 0.2e1 * (t100 * t82 - t101 * t81) * MDP(30) + 0.2e1 * (MDP(18) + MDP(20)) * (-t100 * t92 + t101 * t90); (-MDP(10) * t139 - MDP(9) * t136) * pkin(7) + (-t150 * t137 + t82) * MDP(29) + t141 + t139 * MDP(7) + t136 * MDP(6) + (t101 * t108 - t169) * MDP(20) + (t106 * t92 + t108 * t90) * MDP(23) + (-t100 * t110 - t101 * t109) * MDP(18) + (-t109 * t90 + t110 * t92) * MDP(19) + (t150 * MDP(30) - t158) * t134; MDP(8) + (t109 ^ 2 + t110 ^ 2) * MDP(19) + (t172 + t164) * t108 + 0.2e1 * (-MDP(17) * t135 + t157) * pkin(2) + (MDP(23) * t106 + t144) * t106 + t154; (-t149 * t137 + t82) * MDP(29) + t141 + ((-t100 * t132 - t101 * t133) * MDP(18) + (t132 * t92 - t133 * t90) * MDP(19)) * pkin(3) + (t101 * t124 - t168) * MDP(20) + (t122 * t92 + t124 * t90) * MDP(23) + (t149 * MDP(30) - t158) * t134; (-0.2e1 * pkin(4) - t109) * MDP(21) + (0.2e1 * qJ(5) + t167) * MDP(22) + (t106 * t122 + t108 * t124) * MDP(23) + ((t109 * t133 + t110 * t132) * MDP(19) - t133 * MDP(21) + t132 * MDP(22)) * pkin(3) + (t157 + (MDP(22) * t133 - MDP(17)) * t135) * pkin(2) + t154 + t176 * (t106 + t122); (t132 ^ 2 + t133 ^ 2) * MDP(19) * pkin(3) ^ 2 + (t172 + t163) * t124 + (MDP(23) * t122 + t144) * t122 + t154; MDP(19) * t104 + MDP(23) * t89 - t177 * t101 - t159; 0; 0; MDP(19) + MDP(23); MDP(23) * t90 + (MDP(20) + t147) * t101; MDP(21) + t164; MDP(21) + t163; 0; MDP(23); t80 * MDP(29) - t81 * MDP(30) + t143 + t162; t147 * t105 + t152; t147 * t121 + t152; -t176; t147; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
