% Calculate joint inertia matrix for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR7_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:03
% EndTime: 2019-03-09 05:21:04
% DurationCPUTime: 0.45s
% Computational Cost: add. (777->117), mult. (1303->173), div. (0->0), fcn. (1483->8), ass. (0->62)
t136 = sin(qJ(4));
t137 = sin(qJ(3));
t139 = cos(qJ(4));
t140 = cos(qJ(3));
t116 = t136 * t140 + t139 * t137;
t117 = -t136 * t137 + t139 * t140;
t133 = sin(pkin(10));
t134 = cos(pkin(10));
t98 = t133 * t116 - t134 * t117;
t183 = t98 ^ 2;
t135 = sin(qJ(6));
t138 = cos(qJ(6));
t166 = t135 * MDP(25) + t138 * MDP(26);
t152 = t138 * MDP(28) - t135 * MDP(29);
t153 = t134 * t116 + t133 * t117;
t182 = t153 ^ 2 + t183;
t141 = -pkin(1) - pkin(7);
t171 = pkin(8) - t141;
t157 = t171 * t140;
t158 = t171 * t137;
t144 = t136 * t158 - t139 * t157;
t143 = -t117 * qJ(5) + t144;
t145 = t136 * t157 + t139 * t158;
t92 = -t116 * qJ(5) - t145;
t86 = t133 * t92 - t134 * t143;
t88 = t133 * t143 + t134 * t92;
t181 = t153 * t88 + t86 * t98;
t126 = t139 * pkin(3) + pkin(4);
t172 = pkin(3) * t136;
t108 = t134 * t126 - t133 * t172;
t109 = t133 * t126 + t134 * t172;
t180 = t108 * t98 - t109 * t153;
t179 = -t133 * t153 + t134 * t98;
t151 = -MDP(28) * t135 - MDP(29) * t138;
t123 = t137 * pkin(3) + qJ(2);
t174 = 0.2e1 * t123;
t173 = (pkin(1) * MDP(6));
t170 = MDP(22) * pkin(4);
t85 = t86 * t135;
t169 = t86 * t138;
t168 = t135 * t138;
t163 = t153 * MDP(27);
t160 = MDP(24) * t168;
t131 = t135 ^ 2;
t159 = t131 * MDP(23) + MDP(18) + 0.2e1 * t160;
t105 = t116 * pkin(4) + t123;
t106 = -pkin(5) - t108;
t107 = pkin(9) + t109;
t155 = -t106 * t98 - t107 * t153;
t121 = t133 * pkin(4) + pkin(9);
t122 = -t134 * pkin(4) - pkin(5);
t154 = -t121 * t153 - t122 * t98;
t150 = t117 * MDP(19) - t116 * MDP(20) - t152 * t98;
t149 = (t139 * MDP(19) - t136 * MDP(20)) * pkin(3);
t148 = (MDP(25) * t138 - MDP(26) * t135) * t98;
t132 = t138 ^ 2;
t147 = t117 * MDP(16) - t116 * MDP(17) + t144 * MDP(19) + t145 * MDP(20) + t166 * t153 + (-(-t131 + t132) * MDP(24) - MDP(23) * t168) * t98;
t146 = 0.2e1 * t152;
t89 = pkin(5) * t153 + pkin(9) * t98 + t105;
t84 = t135 * t89 + t138 * t88;
t83 = -t135 * t88 + t138 * t89;
t1 = [MDP(1) + t116 * MDP(19) * t174 + (t105 ^ 2 + t86 ^ 2 + t88 ^ 2) * MDP(22) + (t132 * MDP(23) - 0.2e1 * t160) * t183 + (MDP(7) * t140 - 0.2e1 * t137 * MDP(8)) * t140 + ((-2 * MDP(4) + t173) * pkin(1)) + (MDP(14) * t117 - 0.2e1 * t116 * MDP(15) + MDP(20) * t174) * t117 + (-0.2e1 * t148 + t163) * t153 + (0.2e1 * t137 * MDP(12) + 0.2e1 * t140 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) - 0.2e1 * t181 * MDP(21) + 0.2e1 * (t153 * t83 - t85 * t98) * MDP(28) + 0.2e1 * (-t153 * t84 - t169 * t98) * MDP(29); MDP(4) - t173 + t181 * MDP(22) - (MDP(21) - t151) * t182; t182 * MDP(22) + MDP(6); t180 * MDP(21) + (-t108 * t86 + t109 * t88) * MDP(22) + (t155 * t135 - t169) * MDP(28) + (t155 * t138 + t85) * MDP(29) + (t141 * MDP(12) + MDP(9)) * t140 + (-t141 * MDP(13) - MDP(10)) * t137 + t147; t140 * MDP(12) - t137 * MDP(13) - t180 * MDP(22) + t150; MDP(11) + (t108 ^ 2 + t109 ^ 2) * MDP(22) - t106 * t146 + 0.2e1 * t149 + t159; (t154 * t135 - t169) * MDP(28) + (t154 * t138 + t85) * MDP(29) + (t179 * MDP(21) + (t133 * t88 - t134 * t86) * MDP(22)) * pkin(4) + t147; -t179 * t170 + t150; (t108 * t134 + t109 * t133) * t170 + t149 + t159 - t152 * (t106 + t122); (t133 ^ 2 + t134 ^ 2) * MDP(22) * pkin(4) ^ 2 - t122 * t146 + t159; t105 * MDP(22) + t152 * t153; 0; 0; 0; MDP(22); t83 * MDP(28) - t84 * MDP(29) - t148 + t163; t151 * t153; t151 * t107 + t166; t151 * t121 + t166; t152; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
