% Calculate joint inertia matrix for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR8_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:21
% EndTime: 2019-03-09 02:36:22
% DurationCPUTime: 0.42s
% Computational Cost: add. (602->121), mult. (1082->174), div. (0->0), fcn. (1208->8), ass. (0->69)
t132 = sin(qJ(5));
t135 = cos(qJ(5));
t141 = t132 * MDP(23) + t135 * MDP(24);
t128 = sin(pkin(10));
t129 = cos(pkin(10));
t133 = sin(qJ(4));
t162 = cos(qJ(4));
t108 = t133 * t128 - t162 * t129;
t131 = sin(qJ(6));
t134 = cos(qJ(6));
t112 = t131 * t135 + t134 * t132;
t90 = t112 * t108;
t111 = t131 * t132 - t134 * t135;
t92 = t111 * t108;
t169 = t92 * MDP(27) + t90 * MDP(28);
t163 = pkin(8) + pkin(9);
t115 = t163 * t132;
t116 = t163 * t135;
t145 = t112 * MDP(27) - t111 * MDP(28) + (-t134 * t115 - t131 * t116) * MDP(30) - (-t131 * t115 + t134 * t116) * MDP(31);
t167 = t132 * MDP(20) + t135 * MDP(21) - t141 * pkin(8) + t145;
t166 = t108 ^ 2;
t109 = t162 * t128 + t133 * t129;
t105 = t109 ^ 2;
t165 = -2 * MDP(26);
t164 = 0.2e1 * MDP(31);
t161 = t109 * pkin(5);
t130 = -pkin(1) - qJ(3);
t160 = -pkin(7) + t130;
t89 = t112 * t109;
t91 = t111 * t109;
t159 = -t89 * MDP(30) + t91 * MDP(31);
t113 = t160 * t128;
t114 = t160 * t129;
t95 = t162 * t113 + t133 * t114;
t157 = t135 * t95;
t118 = t128 * pkin(3) + qJ(2);
t93 = t109 * pkin(4) + t108 * pkin(8) + t118;
t80 = t157 + (pkin(9) * t108 + t93) * t132;
t158 = t134 * t80;
t156 = t108 * t132;
t155 = t108 * t135;
t154 = t132 * t135;
t101 = t111 * MDP(30);
t153 = -t112 * MDP(31) - t101;
t117 = t128 ^ 2 + t129 ^ 2;
t152 = MDP(25) * t112;
t151 = t108 * MDP(17);
t148 = MDP(22) + MDP(29);
t147 = t109 * MDP(29) + t169;
t146 = MDP(19) * t154;
t81 = -t132 * t95 + t135 * t93;
t79 = pkin(9) * t155 + t161 + t81;
t76 = -t131 * t80 + t134 * t79;
t144 = t128 * MDP(7) + t129 * MDP(8);
t143 = -t135 * MDP(20) + t132 * MDP(21);
t142 = t135 * MDP(23) - t132 * MDP(24);
t94 = t133 * t113 - t162 * t114;
t140 = -MDP(16) - t142;
t139 = (MDP(30) * t134 - MDP(31) * t131) * pkin(5);
t138 = -t140 + t153;
t136 = qJ(2) ^ 2;
t127 = t135 ^ 2;
t126 = t132 ^ 2;
t121 = -t135 * pkin(5) - pkin(4);
t107 = t117 * t130;
t83 = -pkin(5) * t156 + t94;
t82 = t132 * t93 + t157;
t77 = t131 * t79 + t158;
t1 = [(t117 * t130 ^ 2 + t136) * MDP(10) + (pkin(1) ^ 2 + t136) * MDP(6) - 0.2e1 * pkin(1) * MDP(4) - 0.2e1 * t118 * t151 + MDP(1) + (t92 * MDP(25) - t90 * t165) * t92 + t148 * t105 + (t127 * MDP(18) + MDP(11) - 0.2e1 * t146) * t166 - 0.2e1 * t107 * MDP(9) + 0.2e1 * (t76 * t109 - t83 * t90) * MDP(30) + (-t77 * t109 + t83 * t92) * t164 + 0.2e1 * (t81 * t109 - t94 * t156) * MDP(23) + 0.2e1 * (-t82 * t109 - t94 * t155) * MDP(24) + 0.2e1 * (MDP(5) + t144) * qJ(2) + 0.2e1 * (t118 * MDP(16) + t169 + (MDP(12) + t143) * t108) * t109; MDP(4) - pkin(1) * MDP(6) - t117 * MDP(9) + t107 * MDP(10) + (-t108 * t90 - t89 * t109) * MDP(30) + (t108 * t92 + t91 * t109) * MDP(31) + t141 * (-t105 - t166); t117 * MDP(10) + MDP(6); qJ(2) * MDP(10) + t138 * t109 + t144 - t151; 0; MDP(10); -t95 * MDP(17) + t92 * t152 + (-t92 * t111 + t112 * t90) * MDP(26) + (t83 * t111 - t121 * t90) * MDP(30) + (t83 * t112 + t121 * t92) * MDP(31) + t140 * t94 + (-MDP(13) - MDP(18) * t154 + (t126 - t127) * MDP(19) + t141 * pkin(4)) * t108 + (-MDP(14) + t167) * t109; -t109 * MDP(17) - t138 * t108; 0; 0.2e1 * t146 + 0.2e1 * t121 * t101 + t126 * MDP(18) + MDP(15) + 0.2e1 * t142 * pkin(4) + (t111 * t165 + t121 * t164 + t152) * t112; t109 * MDP(22) + t81 * MDP(23) - t82 * MDP(24) + (t134 * t161 + t76) * MDP(30) + (-t158 + (-t79 - t161) * t131) * MDP(31) + t143 * t108 + t147; -t141 * t109 + t159; t142 + t153; t167; 0.2e1 * t139 + t148; MDP(30) * t76 - MDP(31) * t77 + t147; t159; t153; t145; MDP(29) + t139; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
