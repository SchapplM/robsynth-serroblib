% Calculate joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:26
% EndTime: 2021-01-16 03:28:29
% DurationCPUTime: 0.53s
% Computational Cost: add. (715->132), mult. (1464->204), div. (0->0), fcn. (1726->12), ass. (0->75)
t130 = sin(qJ(6));
t134 = cos(qJ(6));
t143 = t134 * MDP(28) - t130 * MDP(29);
t170 = MDP(21) + t143;
t159 = t130 * MDP(25) + t134 * MDP(26);
t126 = sin(pkin(12));
t128 = cos(pkin(12));
t132 = sin(qJ(3));
t136 = cos(qJ(3));
t113 = t126 * t132 - t128 * t136;
t120 = -pkin(3) * t136 - pkin(2);
t105 = pkin(4) * t113 + t120;
t169 = 0.2e1 * t105;
t168 = pkin(3) * t126;
t166 = qJ(4) + pkin(8);
t131 = sin(qJ(5));
t135 = cos(qJ(5));
t116 = t166 * t132;
t117 = t166 * t136;
t103 = -t128 * t116 - t117 * t126;
t114 = t126 * t136 + t128 * t132;
t91 = -pkin(9) * t114 + t103;
t104 = -t116 * t126 + t117 * t128;
t92 = -pkin(9) * t113 + t104;
t84 = t131 * t92 - t135 * t91;
t165 = t84 * t134;
t102 = -t113 * t131 + t114 * t135;
t164 = t102 * t130;
t163 = t102 * t134;
t127 = sin(pkin(6));
t133 = sin(qJ(2));
t162 = t127 * t133;
t137 = cos(qJ(2));
t161 = t127 * t137;
t160 = t130 * t134;
t158 = MDP(10) * t136;
t101 = t113 * t135 + t114 * t131;
t157 = t101 * MDP(27);
t156 = t102 * MDP(22);
t119 = pkin(3) * t128 + pkin(4);
t108 = t119 * t135 - t131 * t168;
t155 = t108 * MDP(21);
t109 = -t119 * t131 - t135 * t168;
t154 = t109 * MDP(22);
t153 = t113 * MDP(12);
t152 = t114 * MDP(13);
t151 = t120 * MDP(15);
t148 = MDP(24) * t160;
t124 = t130 ^ 2;
t147 = t124 * MDP(23) + MDP(20) + 0.2e1 * t148;
t146 = -pkin(5) * t102 - pkin(10) * t101;
t106 = -pkin(5) - t108;
t107 = pkin(10) - t109;
t145 = -t101 * t107 + t102 * t106;
t144 = MDP(25) * t134 - MDP(26) * t130;
t142 = -MDP(28) * t130 - MDP(29) * t134;
t129 = cos(pkin(6));
t111 = t129 * t136 - t132 * t162;
t112 = t129 * t132 + t136 * t162;
t94 = t111 * t128 - t112 * t126;
t95 = t111 * t126 + t112 * t128;
t88 = t131 * t95 - t135 * t94;
t89 = t131 * t94 + t135 * t95;
t141 = -t89 * MDP(22) - t170 * t88;
t125 = t134 ^ 2;
t85 = t131 * t91 + t135 * t92;
t140 = -t84 * MDP(21) - t85 * MDP(22) + ((-t124 + t125) * MDP(24) + MDP(23) * t160 + MDP(18)) * t102 + (-MDP(19) + t159) * t101;
t139 = t151 + t152 + t153 + t156;
t83 = pkin(5) * t101 - pkin(10) * t102 + t105;
t79 = t84 * t130;
t78 = -t130 * t161 + t134 * t89;
t77 = -t130 * t89 - t134 * t161;
t76 = t130 * t83 + t134 * t85;
t75 = -t130 * t85 + t134 * t83;
t1 = [MDP(1) + (t127 ^ 2 * t137 ^ 2 + t94 ^ 2 + t95 ^ 2) * MDP(15); (-t113 * t95 - t114 * t94) * MDP(14) + (t103 * t94 + t104 * t95) * MDP(15) + (t101 * t77 + t164 * t88) * MDP(28) + (-t101 * t78 + t163 * t88) * MDP(29) + (-t133 * MDP(4) + (-MDP(11) * t132 - MDP(21) * t101 + MDP(3) - t139 + t158) * t137) * t127; MDP(2) + 0.2e1 * pkin(2) * t158 + (t103 ^ 2 + t104 ^ 2) * MDP(15) + t156 * t169 + (t151 + 0.2e1 * t152 + 0.2e1 * t153) * t120 + (-0.2e1 * MDP(11) * pkin(2) + MDP(5) * t132 + 0.2e1 * MDP(6) * t136) * t132 + (MDP(23) * t125 + MDP(16) - 0.2e1 * t148) * t102 ^ 2 + (MDP(21) * t169 + t157 + 0.2e1 * (-MDP(17) + t144) * t102) * t101 + 0.2e1 * (-t103 * t114 - t104 * t113) * MDP(14) + 0.2e1 * (t101 * t75 + t164 * t84) * MDP(28) + 0.2e1 * (-t101 * t76 + t163 * t84) * MDP(29); t111 * MDP(10) - t112 * MDP(11) + t94 * MDP(12) - t95 * MDP(13) + (t126 * t95 + t128 * t94) * MDP(15) * pkin(3) + t141; (-MDP(10) * t132 - MDP(11) * t136) * pkin(8) + ((-t113 * t126 - t114 * t128) * MDP(14) + (t103 * t128 + t104 * t126) * MDP(15)) * pkin(3) + t140 + t103 * MDP(12) - t104 * MDP(13) + t132 * MDP(7) + t136 * MDP(8) + (t130 * t145 - t165) * MDP(28) + (t134 * t145 + t79) * MDP(29); MDP(9) + (t126 ^ 2 + t128 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * t155 + 0.2e1 * t154 + t147 - 0.2e1 * t143 * t106 + 0.2e1 * (t128 * MDP(12) - t126 * MDP(13)) * pkin(3); -MDP(15) * t161; t101 * t170 + t139; 0; MDP(15); t141; (t130 * t146 - t165) * MDP(28) + (t134 * t146 + t79) * MDP(29) + t140; t147 + t154 + t155 + t143 * (pkin(5) - t106); 0; 0.2e1 * pkin(5) * t143 + t147; MDP(28) * t77 - MDP(29) * t78; t75 * MDP(28) - t76 * MDP(29) + t102 * t144 + t157; t107 * t142 + t159; t143; pkin(10) * t142 + t159; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
