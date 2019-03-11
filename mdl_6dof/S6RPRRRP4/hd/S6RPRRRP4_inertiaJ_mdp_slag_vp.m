% Calculate joint inertia matrix for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP4_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:47
% EndTime: 2019-03-09 06:08:48
% DurationCPUTime: 0.47s
% Computational Cost: add. (1097->148), mult. (2015->202), div. (0->0), fcn. (2331->8), ass. (0->73)
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t163 = t135 * MDP(24) + t138 * MDP(25);
t157 = MDP(28) * t135;
t133 = sin(pkin(10));
t134 = cos(pkin(10));
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t112 = t133 * t137 - t134 * t140;
t121 = -pkin(2) * t134 - pkin(1);
t106 = pkin(3) * t112 + t121;
t176 = 0.2e1 * t106;
t175 = 0.2e1 * t121;
t174 = 2 * MDP(29);
t173 = pkin(1) * MDP(7);
t139 = cos(qJ(4));
t172 = pkin(3) * t139;
t123 = -pkin(4) - t172;
t171 = pkin(4) - t123;
t170 = pkin(7) + qJ(2);
t169 = MDP(30) * pkin(5);
t136 = sin(qJ(4));
t113 = t133 * t140 + t134 * t137;
t115 = t170 * t133;
t116 = t170 * t134;
t151 = -t140 * t115 - t116 * t137;
t96 = -pkin(8) * t113 + t151;
t146 = t115 * t137 - t116 * t140;
t97 = -pkin(8) * t112 - t146;
t93 = t136 * t96 + t139 * t97;
t168 = t138 * t93;
t167 = t133 * MDP(5);
t166 = t134 * MDP(4);
t165 = t135 * t138;
t128 = t138 * qJ(6);
t92 = t136 * t97 - t139 * t96;
t164 = t92 * MDP(27);
t131 = t135 ^ 2;
t132 = t138 ^ 2;
t161 = t131 + t132;
t160 = MDP(25) * t135;
t104 = t139 * t112 + t113 * t136;
t159 = MDP(26) * t104;
t158 = MDP(27) * t138;
t105 = -t112 * t136 + t113 * t139;
t156 = MDP(29) * t105;
t155 = t112 * MDP(13);
t154 = t135 * MDP(29);
t124 = -pkin(5) * t138 - pkin(4);
t153 = MDP(23) * t165;
t152 = t131 * MDP(22) + MDP(19) + 0.2e1 * t153;
t94 = pkin(4) * t104 - pkin(9) * t105 + t106;
t85 = -t135 * t93 + t138 * t94;
t150 = -pkin(4) * t105 - pkin(9) * t104;
t82 = pkin(5) * t104 - t105 * t128 + t85;
t84 = t168 + (-qJ(6) * t105 + t94) * t135;
t149 = t135 * t84 + t138 * t82;
t148 = t85 * MDP(27) - (t135 * t94 + t168) * MDP(28);
t122 = pkin(3) * t136 + pkin(9);
t147 = -t104 * t122 + t105 * t123;
t145 = -t157 + t158;
t144 = t135 * MDP(27) + t138 * MDP(28);
t143 = (MDP(20) * t139 - MDP(21) * t136) * pkin(3);
t142 = t84 * t138 * MDP(29) - t93 * MDP(21) + (t157 - MDP(20)) * t92 + (MDP(17) + (-t131 + t132) * MDP(23) + MDP(22) * t165) * t105 + (-MDP(18) + t163) * t104;
t119 = pkin(9) * t138 + t128;
t118 = (-qJ(6) - pkin(9)) * t135;
t117 = t124 - t172;
t114 = t119 * t138;
t110 = t122 * t138 + t128;
t109 = (-qJ(6) - t122) * t135;
t107 = t110 * t138;
t87 = pkin(5) * t105 * t135 + t92;
t1 = [MDP(1) + t155 * t175 + (t82 ^ 2 + t84 ^ 2 + t87 ^ 2) * MDP(30) + (0.2e1 * t166 - 0.2e1 * t167 + t173) * pkin(1) + (MDP(14) * t175 + MDP(8) * t113 - 0.2e1 * t112 * MDP(9)) * t113 + (MDP(21) * t176 - 0.2e1 * t149 * MDP(29) + 0.2e1 * t144 * t92 + (t132 * MDP(22) + MDP(15) - 0.2e1 * t153) * t105) * t105 + (MDP(20) * t176 + t159 + 0.2e1 * (MDP(24) * t138 - MDP(16) - t160) * t105 + 0.2e1 * t148) * t104 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t133 ^ 2 + t134 ^ 2) * qJ(2); -t166 + t167 - t173 + t155 + t113 * MDP(14) + t149 * MDP(30) + (-t161 * MDP(29) + MDP(21)) * t105 + (MDP(20) + t145) * t104; t161 * MDP(30) + MDP(7); t146 * MDP(14) + t151 * MDP(13) - t112 * MDP(11) + t113 * MDP(10) + (t109 * t82 + t110 * t84 + t117 * t87) * MDP(30) + t142 + (t147 * MDP(27) + (-t105 * t110 - t82) * MDP(29)) * t135 + (t147 * MDP(28) - t109 * t156 - t164) * t138; (t109 * t138 + t110 * t135) * MDP(30); MDP(12) + (-t109 * t135 + t107) * t174 + (t109 ^ 2 + t110 ^ 2 + t117 ^ 2) * MDP(30) + t152 - 0.2e1 * t145 * t123 + 0.2e1 * t143; (t118 * t82 + t119 * t84 + t124 * t87) * MDP(30) + (t150 * MDP(28) - t118 * t156 - t164) * t138 + (t150 * MDP(27) + (-t105 * t119 - t82) * MDP(29)) * t135 + t142; (t118 * t138 + t119 * t135) * MDP(30); (t107 + t114) * MDP(29) + (t109 * t118 + t110 * t119 + t117 * t124) * MDP(30) + t171 * t158 + t143 + (-t171 * MDP(28) + (-t109 - t118) * MDP(29)) * t135 + t152; (-t118 * t135 + t114) * t174 + (t118 ^ 2 + t119 ^ 2 + t124 ^ 2) * MDP(30) + 0.2e1 * t145 * pkin(4) + t152; t82 * t169 + t159 + (-t160 + (-MDP(29) * pkin(5) + MDP(24)) * t138) * t105 + t148; -t157 + (MDP(27) + t169) * t138; -t144 * t122 + (t109 * MDP(30) - t154) * pkin(5) + t163; -t144 * pkin(9) + (t118 * MDP(30) - t154) * pkin(5) + t163; MDP(30) * pkin(5) ^ 2 + MDP(26); t87 * MDP(30); 0; t117 * MDP(30); t124 * MDP(30); 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
