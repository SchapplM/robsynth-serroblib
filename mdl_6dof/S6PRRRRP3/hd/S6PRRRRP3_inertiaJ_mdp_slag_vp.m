% Calculate joint inertia matrix for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP3_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:09
% EndTime: 2019-03-09 00:11:11
% DurationCPUTime: 0.62s
% Computational Cost: add. (682->166), mult. (1443->250), div. (0->0), fcn. (1538->10), ass. (0->75)
t124 = sin(qJ(4));
t128 = cos(qJ(4));
t123 = sin(qJ(5));
t127 = cos(qJ(5));
t146 = t127 * t128;
t104 = t123 * t124 - t146;
t105 = t123 * t128 + t124 * t127;
t160 = pkin(9) + pkin(10);
t109 = t160 * t124;
t110 = t160 * t128;
t90 = -t127 * t109 - t110 * t123;
t91 = -t109 * t123 + t110 * t127;
t137 = t105 * MDP(21) - t104 * MDP(22) + t90 * MDP(24) - t91 * MDP(25);
t166 = t137 - (MDP(17) * t124 + MDP(18) * t128) * pkin(9) + t124 * MDP(14) + t128 * MDP(15);
t141 = t127 * MDP(24);
t165 = (-MDP(25) * t123 + t141) * pkin(4);
t125 = sin(qJ(3));
t97 = t105 * t125;
t149 = t124 * t125;
t98 = -t123 * t149 + t125 * t146;
t154 = t98 * MDP(21) - t97 * MDP(22);
t129 = cos(qJ(3));
t108 = -pkin(3) * t129 - pkin(9) * t125 - pkin(2);
t103 = t128 * t108;
t147 = t125 * t128;
t158 = pkin(8) * t124;
t83 = -pkin(10) * t147 + t103 + (-pkin(4) - t158) * t129;
t156 = pkin(8) * t129;
t139 = t128 * t156;
t87 = t139 + (-pkin(10) * t125 + t108) * t124;
t76 = -t123 * t87 + t127 * t83;
t77 = t123 * t83 + t127 * t87;
t164 = t76 * MDP(24) - t77 * MDP(25) + t154;
t162 = -2 * MDP(20);
t161 = 2 * MDP(26);
t159 = pkin(4) * t123;
t157 = pkin(8) * t128;
t122 = cos(pkin(6));
t121 = sin(pkin(6));
t126 = sin(qJ(2));
t151 = t121 * t126;
t100 = t122 * t125 + t129 * t151;
t130 = cos(qJ(2));
t150 = t121 * t130;
t85 = -t100 * t124 - t128 * t150;
t86 = t100 * t128 - t124 * t150;
t78 = -t123 * t86 + t127 * t85;
t79 = t123 * t85 + t127 * t86;
t155 = t78 * MDP(24) - t79 * MDP(25);
t153 = MDP(27) * pkin(5);
t152 = MDP(19) * t98;
t148 = t124 * t128;
t115 = -pkin(4) * t128 - pkin(3);
t94 = pkin(5) * t104 + t115;
t143 = t94 * MDP(27);
t107 = pkin(4) * t149 + t125 * pkin(8);
t142 = MDP(11) * t125;
t140 = MDP(16) + MDP(23);
t138 = MDP(13) * t148;
t136 = MDP(14) * t128 - MDP(15) * t124;
t134 = MDP(17) * t128 - MDP(18) * t124;
t132 = MDP(24) * t104 + MDP(25) * t105;
t119 = t128 ^ 2;
t118 = t125 ^ 2;
t117 = t124 ^ 2;
t114 = pkin(4) * t127 + pkin(5);
t99 = -t122 * t129 + t125 * t151;
t93 = t108 * t124 + t139;
t92 = -t124 * t156 + t103;
t84 = pkin(5) * t97 + t107;
t81 = -qJ(6) * t104 + t91;
t80 = -qJ(6) * t105 + t90;
t73 = -qJ(6) * t97 + t77;
t72 = -pkin(5) * t129 - qJ(6) * t98 + t76;
t1 = [MDP(1) + (t78 ^ 2 + t79 ^ 2 + t99 ^ 2) * MDP(27); (-t129 * t85 + t99 * t149) * MDP(17) + (t129 * t86 + t99 * t147) * MDP(18) + (-t129 * t78 + t97 * t99) * MDP(24) + (t129 * t79 + t98 * t99) * MDP(25) + (-t78 * t98 - t79 * t97) * MDP(26) + (t72 * t78 + t73 * t79 + t84 * t99) * MDP(27) + (-t126 * MDP(4) + (MDP(10) * t129 + MDP(3) - t142) * t130) * t121; MDP(2) - 0.2e1 * pkin(2) * t142 + (t72 ^ 2 + t73 ^ 2 + t84 ^ 2) * MDP(27) + (t97 * t162 + t152) * t98 + t140 * t129 ^ 2 + (MDP(12) * t119 + MDP(5) - 0.2e1 * t138) * t118 + 0.2e1 * (MDP(10) * pkin(2) + (MDP(6) - t136) * t125 - t154) * t129 + 0.2e1 * (t118 * t158 - t129 * t92) * MDP(17) + 0.2e1 * (t118 * t157 + t129 * t93) * MDP(18) + 0.2e1 * (t107 * t97 - t129 * t76) * MDP(24) + 0.2e1 * (t107 * t98 + t129 * t77) * MDP(25) + (-t72 * t98 - t73 * t97) * t161; -t100 * MDP(11) + (-t104 * t79 - t105 * t78) * MDP(26) + (t78 * t80 + t79 * t81) * MDP(27) + (-MDP(10) + t132 - t134 + t143) * t99; t105 * t152 + (-t104 * t98 - t105 * t97) * MDP(20) + (t104 * t107 + t115 * t97) * MDP(24) + (t105 * t107 + t115 * t98) * MDP(25) + (-t104 * t73 - t105 * t72 - t80 * t98 - t81 * t97) * MDP(26) + (t72 * t80 + t73 * t81 + t84 * t94) * MDP(27) + (-pkin(8) * MDP(11) + MDP(8) - t166) * t129 + (MDP(7) - pkin(8) * MDP(10) + MDP(12) * t148 + (-t117 + t119) * MDP(13) + (-pkin(3) * t124 - t157) * MDP(17) + (-pkin(3) * t128 + t158) * MDP(18)) * t125; MDP(9) + t117 * MDP(12) + 0.2e1 * t138 - t104 * t81 * t161 + (t80 ^ 2 + t81 ^ 2 + t94 ^ 2) * MDP(27) + (MDP(19) * t105 + t104 * t162 - t80 * t161) * t105 + 0.2e1 * t134 * pkin(3) + 0.2e1 * t132 * t115; t85 * MDP(17) - t86 * MDP(18) + (t114 * t78 + t79 * t159) * MDP(27) + t155; t92 * MDP(17) - t93 * MDP(18) + (-t114 * t98 - t97 * t159) * MDP(26) + (t114 * t72 + t73 * t159) * MDP(27) + (-t140 - t165) * t129 + t136 * t125 + t164; (-t104 * t159 - t105 * t114) * MDP(26) + (t114 * t80 + t81 * t159) * MDP(27) + t166; t114 ^ 2 * MDP(27) + (0.2e1 * t141 + (MDP(27) * t159 - 0.2e1 * MDP(25)) * t123) * pkin(4) + t140; t78 * t153 + t155; -t129 * MDP(23) + (-MDP(26) * t98 + MDP(27) * t72) * pkin(5) + t164; (-MDP(26) * t105 + MDP(27) * t80) * pkin(5) + t137; t114 * t153 + MDP(23) + t165; MDP(27) * pkin(5) ^ 2 + MDP(23); t99 * MDP(27); t84 * MDP(27); t143; 0; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
