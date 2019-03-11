% Calculate joint inertia matrix for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:31
% EndTime: 2019-03-08 19:01:32
% DurationCPUTime: 0.38s
% Computational Cost: add. (454->111), mult. (1094->179), div. (0->0), fcn. (1297->14), ass. (0->71)
t134 = sin(qJ(6));
t138 = cos(qJ(6));
t159 = t134 * MDP(22) + t138 * MDP(23);
t148 = t138 * MDP(25) - t134 * MDP(26);
t170 = -MDP(18) - t148;
t140 = cos(qJ(4));
t121 = -t140 * pkin(4) - pkin(3);
t169 = 0.2e1 * t121;
t168 = pkin(9) + pkin(10);
t136 = sin(qJ(4));
t115 = t168 * t136;
t116 = t168 * t140;
t135 = sin(qJ(5));
t139 = cos(qJ(5));
t102 = t139 * t115 + t135 * t116;
t166 = t102 * t138;
t114 = t135 * t140 + t139 * t136;
t165 = t114 * t134;
t164 = t114 * t138;
t129 = sin(pkin(7));
t137 = sin(qJ(3));
t163 = t129 * t137;
t141 = cos(qJ(3));
t162 = t129 * t141;
t131 = cos(pkin(13));
t132 = cos(pkin(7));
t161 = t131 * t132;
t160 = t134 * t138;
t113 = t135 * t136 - t139 * t140;
t158 = t113 * MDP(24);
t157 = t114 * MDP(19);
t154 = t140 * MDP(11);
t153 = MDP(21) * t160;
t126 = t134 ^ 2;
t152 = t126 * MDP(20) + MDP(17) + 0.2e1 * t153;
t151 = -pkin(5) * t114 - pkin(11) * t113;
t119 = t135 * pkin(4) + pkin(11);
t120 = -t139 * pkin(4) - pkin(5);
t150 = -t113 * t119 + t114 * t120;
t149 = MDP(22) * t138 - MDP(23) * t134;
t147 = -MDP(25) * t134 - MDP(26) * t138;
t130 = sin(pkin(6));
t133 = cos(pkin(6));
t107 = -t130 * t131 * t129 + t133 * t132;
t128 = sin(pkin(13));
t98 = t133 * t163 + (t128 * t141 + t137 * t161) * t130;
t86 = t107 * t140 - t98 * t136;
t87 = t107 * t136 + t98 * t140;
t82 = t135 * t87 - t139 * t86;
t83 = t135 * t86 + t139 * t87;
t146 = -t83 * MDP(19) + t170 * t82;
t108 = t140 * t132 - t136 * t163;
t109 = t136 * t132 + t140 * t163;
t93 = -t139 * t108 + t135 * t109;
t94 = t135 * t108 + t139 * t109;
t145 = -t94 * MDP(19) + t170 * t93;
t144 = (t139 * MDP(18) - t135 * MDP(19)) * pkin(4);
t103 = -t135 * t115 + t139 * t116;
t127 = t138 ^ 2;
t143 = -t102 * MDP(18) - t103 * MDP(19) + (MDP(20) * t160 + MDP(15) + (-t126 + t127) * MDP(21)) * t114 + (-MDP(16) + t159) * t113;
t142 = -t136 * MDP(12) - t113 * MDP(18) + MDP(4) + t154 - t157;
t99 = t102 * t134;
t97 = -t133 * t162 + (t128 * t137 - t141 * t161) * t130;
t96 = t113 * pkin(5) - t114 * pkin(11) + t121;
t89 = -t134 * t162 + t138 * t94;
t88 = -t134 * t94 - t138 * t162;
t85 = t138 * t103 + t134 * t96;
t84 = -t134 * t103 + t138 * t96;
t78 = t97 * t134 + t138 * t83;
t77 = -t134 * t83 + t97 * t138;
t1 = [MDP(1) + (t133 ^ 2 + (t128 ^ 2 + t131 ^ 2) * t130 ^ 2) * MDP(2); t133 * MDP(2); MDP(2); -t98 * MDP(5) + (t77 * t113 + t82 * t165) * MDP(25) + (-t78 * t113 + t82 * t164) * MDP(26) - t142 * t97; (t88 * t113 + t93 * t165) * MDP(25) + (-t89 * t113 + t93 * t164) * MDP(26) + (-t137 * MDP(5) + t142 * t141) * t129; 0.2e1 * pkin(3) * t154 + t157 * t169 + MDP(3) + (MDP(18) * t169 + t158 + 0.2e1 * (-MDP(14) + t149) * t114) * t113 + 0.2e1 * (t102 * t165 + t84 * t113) * MDP(25) + 0.2e1 * (t102 * t164 - t85 * t113) * MDP(26) + (-0.2e1 * pkin(3) * MDP(12) + MDP(6) * t136 + 0.2e1 * t140 * MDP(7)) * t136 + (t127 * MDP(20) + MDP(13) - 0.2e1 * t153) * t114 ^ 2; t86 * MDP(11) - t87 * MDP(12) + t146; t108 * MDP(11) - t109 * MDP(12) + t145; t136 * MDP(8) + t140 * MDP(9) + (t150 * t134 - t166) * MDP(25) + (t150 * t138 + t99) * MDP(26) + (-t136 * MDP(11) - t140 * MDP(12)) * pkin(9) + t143; -0.2e1 * t120 * t148 + MDP(10) + 0.2e1 * t144 + t152; t146; t145; (t151 * t134 - t166) * MDP(25) + (t151 * t138 + t99) * MDP(26) + t143; t144 + t152 + t148 * (pkin(5) - t120); 0.2e1 * pkin(5) * t148 + t152; t77 * MDP(25) - t78 * MDP(26); t88 * MDP(25) - t89 * MDP(26); t84 * MDP(25) - t85 * MDP(26) + t149 * t114 + t158; t147 * t119 + t159; t147 * pkin(11) + t159; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
