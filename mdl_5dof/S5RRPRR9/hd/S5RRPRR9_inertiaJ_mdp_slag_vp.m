% Calculate joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR9_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:59
% EndTime: 2019-12-31 20:22:01
% DurationCPUTime: 0.41s
% Computational Cost: add. (502->107), mult. (964->168), div. (0->0), fcn. (1049->8), ass. (0->60)
t143 = qJ(3) + pkin(6);
t109 = sin(pkin(9));
t110 = cos(pkin(9));
t113 = sin(qJ(2));
t136 = cos(qJ(2));
t96 = t109 * t136 + t110 * t113;
t111 = sin(qJ(5));
t112 = sin(qJ(4));
t114 = cos(qJ(5));
t115 = cos(qJ(4));
t99 = t111 * t115 + t112 * t114;
t75 = t99 * t96;
t98 = t111 * t112 - t114 * t115;
t76 = t98 * t96;
t142 = -t76 * MDP(22) - t75 * MDP(23);
t103 = pkin(2) * t109 + pkin(7);
t120 = MDP(18) * t112 + MDP(19) * t115;
t135 = pkin(8) + t103;
t91 = t135 * t112;
t92 = t135 * t115;
t122 = t99 * MDP(22) - t98 * MDP(23) + (-t111 * t92 - t114 * t91) * MDP(25) - (-t111 * t91 + t114 * t92) * MDP(26);
t141 = t112 * MDP(15) + t115 * MDP(16) - t120 * t103 + t122;
t121 = MDP(18) * t115 - MDP(19) * t112;
t87 = t98 * MDP(25);
t134 = -t99 * MDP(26) - t87;
t140 = t121 + t134;
t139 = -2 * MDP(21);
t138 = 0.2e1 * MDP(26);
t95 = t109 * t113 - t110 * t136;
t137 = pkin(4) * t95;
t133 = t112 * t96;
t101 = t143 * t136;
t123 = t143 * t113;
t85 = t110 * t101 - t109 * t123;
t131 = t115 * t85;
t106 = -t136 * pkin(2) - pkin(1);
t82 = t95 * pkin(3) - t96 * pkin(7) + t106;
t69 = t131 + (-pkin(8) * t96 + t82) * t112;
t132 = t114 * t69;
t130 = t115 * t96;
t129 = MDP(20) * t99;
t128 = t112 * t115;
t127 = 0.2e1 * t136;
t126 = MDP(17) + MDP(24);
t125 = t95 * MDP(24) + t142;
t104 = -pkin(2) * t110 - pkin(3);
t124 = MDP(14) * t128;
t70 = -t112 * t85 + t115 * t82;
t68 = -pkin(8) * t130 + t137 + t70;
t65 = -t111 * t69 + t114 * t68;
t83 = t101 * t109 + t110 * t123;
t119 = (MDP(25) * t114 - MDP(26) * t111) * pkin(4);
t118 = (MDP(15) * t115 - MDP(16) * t112) * t96;
t108 = t115 ^ 2;
t107 = t112 ^ 2;
t100 = -pkin(4) * t115 + t104;
t72 = pkin(4) * t133 + t83;
t71 = t112 * t82 + t131;
t66 = t111 * t68 + t132;
t1 = [MDP(1) + pkin(1) * MDP(9) * t127 + (t106 ^ 2 + t83 ^ 2 + t85 ^ 2) * MDP(12) + (MDP(13) * t108 - 0.2e1 * t124) * t96 ^ 2 + t126 * t95 ^ 2 - (-MDP(20) * t76 + t75 * t139) * t76 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t113 + MDP(5) * t127) * t113 + 0.2e1 * (t118 + t142) * t95 + 0.2e1 * (t83 * t96 - t85 * t95) * MDP(11) + 0.2e1 * (t83 * t133 + t70 * t95) * MDP(18) + 0.2e1 * (t83 * t130 - t71 * t95) * MDP(19) + 0.2e1 * (t65 * t95 + t72 * t75) * MDP(25) + (-t66 * t95 - t72 * t76) * t138; t113 * MDP(6) + t136 * MDP(7) - t76 * t129 + (-t75 * t99 + t76 * t98) * MDP(21) + (t100 * t75 + t72 * t98) * MDP(25) + (-t100 * t76 + t72 * t99) * MDP(26) - t121 * t83 + (-t136 * MDP(10) - t113 * MDP(9)) * pkin(6) + (MDP(13) * t128 + (-t107 + t108) * MDP(14) + t120 * t104) * t96 + t141 * t95 + ((-t109 * t95 - t110 * t96) * MDP(11) + (t109 * t85 - t110 * t83) * MDP(12)) * pkin(2); 0.2e1 * t124 + 0.2e1 * t100 * t87 + t107 * MDP(13) + MDP(8) + (t109 ^ 2 + t110 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t121 * t104 + (t100 * t138 + t98 * t139 + t129) * t99; MDP(12) * t106 + t140 * t95; 0; MDP(12); t95 * MDP(17) + t70 * MDP(18) - t71 * MDP(19) + (t114 * t137 + t65) * MDP(25) + (-t132 + (-t68 - t137) * t111) * MDP(26) + t118 + t125; t141; t140; 0.2e1 * t119 + t126; t65 * MDP(25) - t66 * MDP(26) + t125; t122; t134; MDP(24) + t119; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
