% Calculate joint inertia matrix for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:08
% EndTime: 2019-03-09 08:12:09
% DurationCPUTime: 0.46s
% Computational Cost: add. (764->135), mult. (1345->192), div. (0->0), fcn. (1463->8), ass. (0->69)
t124 = sin(pkin(10));
t126 = cos(pkin(10));
t150 = t124 ^ 2 + t126 ^ 2;
t165 = MDP(20) * t150;
t164 = t150 * MDP(19) - MDP(14);
t128 = sin(qJ(6));
t130 = cos(qJ(6));
t108 = -t124 * t128 + t126 * t130;
t127 = cos(pkin(9));
t120 = -pkin(2) * t127 - pkin(3);
t163 = MDP(16) * t120;
t125 = sin(pkin(9));
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t105 = t125 * t129 - t127 * t131;
t107 = t125 * t131 + t127 * t129;
t121 = -pkin(2) * t131 - pkin(1);
t135 = -qJ(4) * t107 + t121;
t96 = pkin(3) * t105 + t135;
t162 = -0.2e1 * t96;
t117 = pkin(2) * t125 + qJ(4);
t110 = pkin(5) * t124 + t117;
t161 = 0.2e1 * t110;
t160 = 0.2e1 * t131;
t159 = -2 * MDP(22);
t116 = -qJ(5) + t120;
t158 = -pkin(8) + t116;
t157 = -qJ(3) - pkin(7);
t87 = (pkin(3) + qJ(5)) * t105 + t135;
t112 = t157 * t129;
t113 = t157 * t131;
t97 = -t127 * t112 - t113 * t125;
t90 = pkin(4) * t107 + t97;
t84 = t124 * t90 + t126 * t87;
t106 = t124 * t130 + t126 * t128;
t93 = t106 * t105;
t156 = MDP(23) * t93;
t92 = t108 * t105;
t155 = MDP(24) * t92;
t154 = t105 * t117;
t153 = t105 * t126;
t149 = MDP(21) * t108;
t148 = MDP(25) * t107;
t147 = MDP(26) * t106;
t146 = t124 * MDP(17);
t145 = t126 * MDP(18);
t144 = MDP(16) + t165;
t143 = MDP(11) + MDP(13);
t99 = t112 * t125 - t113 * t127;
t142 = t97 ^ 2 + t99 ^ 2;
t89 = t126 * t90;
t83 = -t124 * t87 + t89;
t80 = t124 * t84 + t126 * t83;
t141 = -t124 * t83 + t126 * t84;
t81 = pkin(5) * t107 + t89 + (-pkin(8) * t105 - t87) * t124;
t82 = pkin(8) * t153 + t84;
t140 = (-t128 * t82 + t130 * t81) * MDP(26) - (t128 * t81 + t130 * t82) * MDP(27);
t139 = -MDP(26) * t92 + MDP(27) * t93;
t138 = MDP(17) * t126 - MDP(18) * t124;
t137 = MDP(26) * t108 - t106 * MDP(27);
t136 = -MDP(27) * t108 - t147;
t103 = t158 * t124;
t104 = t158 * t126;
t134 = t108 * MDP(23) - t106 * MDP(24) + (-t103 * t128 + t104 * t130) * MDP(26) - (t103 * t130 + t104 * t128) * MDP(27);
t133 = t136 - t145 - t146;
t101 = t150 * t116;
t91 = -pkin(4) * t105 + t99;
t85 = (-pkin(5) * t126 - pkin(4)) * t105 + t99;
t1 = [MDP(1) + pkin(1) * MDP(9) * t160 + (t121 ^ 2 + t142) * MDP(12) + t105 * MDP(14) * t162 + (t96 ^ 2 + t142) * MDP(16) + (t83 ^ 2 + t84 ^ 2 + t91 ^ 2) * MDP(20) + (MDP(21) * t93 - t92 * t159) * t93 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t129 + MDP(5) * t160) * t129 + (MDP(15) * t162 + t148 + 0.2e1 * t155 + 0.2e1 * t156) * t107 + 0.2e1 * t139 * t85 + 0.2e1 * (t83 * MDP(17) - t84 * MDP(18) + t143 * t97 + t140) * t107 + 0.2e1 * (t141 * MDP(19) - t138 * t91 - t143 * t99) * t105; t129 * MDP(6) + t131 * MDP(7) - MDP(13) * t154 + t97 * MDP(14) + t99 * MDP(15) + (t117 * t99 + t120 * t97) * MDP(16) + (-t117 * t153 + t91 * t124) * MDP(17) + (t124 * t154 + t126 * t91) * MDP(18) - t80 * MDP(19) + (t80 * t116 + t117 * t91) * MDP(20) + t93 * t149 + (-t106 * t93 + t108 * t92) * MDP(22) + (t106 * t85 - t110 * t92) * MDP(26) + (t108 * t85 + t110 * t93) * MDP(27) + (t120 * MDP(13) + t138 * t116 + t134) * t107 + (-MDP(10) * t131 - MDP(9) * t129) * pkin(7) + ((-t105 * t125 - t107 * t127) * MDP(11) + (t125 * t99 - t127 * t97) * MDP(12)) * pkin(2); t147 * t161 + MDP(8) + (t125 ^ 2 + t127 ^ 2) * MDP(12) * pkin(2) ^ 2 + t116 ^ 2 * t165 - 0.2e1 * MDP(19) * t101 + (MDP(27) * t161 + t106 * t159 + t149) * t108 + (0.2e1 * MDP(15) + 0.2e1 * t146 + 0.2e1 * t145 + (MDP(16) + MDP(20)) * t117) * t117 + (0.2e1 * MDP(14) + t163) * t120; t121 * MDP(12) + t96 * MDP(16) + t141 * MDP(20) + t164 * t105 + (-MDP(15) + t133) * t107; 0; MDP(12) + t144; MDP(16) * t97 + MDP(20) * t80 + (MDP(13) + t137 + t138) * t107; MDP(20) * t101 + t163 - t164; 0; t144; MDP(20) * t91 - t138 * t105 + t139; MDP(20) * t117 - t133; 0; 0; MDP(20); t140 + t148 + t155 + t156; t134; t136; t137; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
