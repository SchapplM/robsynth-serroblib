% Calculate joint inertia matrix for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP6_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:09
% EndTime: 2021-01-15 22:37:13
% DurationCPUTime: 0.54s
% Computational Cost: add. (700->161), mult. (1323->230), div. (0->0), fcn. (1254->6), ass. (0->59)
t146 = MDP(21) + MDP(25);
t112 = sin(qJ(3));
t114 = cos(qJ(3));
t145 = -(t112 * MDP(16) + t114 * MDP(17)) * pkin(7) + t112 * MDP(13) + t114 * MDP(14);
t144 = 2 * MDP(18);
t143 = 2 * MDP(19);
t142 = 2 * MDP(20);
t141 = 2 * MDP(22);
t140 = 2 * MDP(23);
t139 = 2 * MDP(24);
t138 = pkin(6) * t112;
t137 = pkin(6) * t114;
t115 = cos(qJ(2));
t136 = pkin(6) * t115;
t135 = qJ(4) + pkin(7);
t110 = sin(pkin(8));
t111 = cos(pkin(8));
t113 = sin(qJ(2));
t132 = t113 * t114;
t97 = -t115 * pkin(2) - t113 * pkin(7) - pkin(1);
t94 = t114 * t97;
t81 = -qJ(4) * t132 + t94 + (-pkin(3) - t138) * t115;
t126 = t114 * t136;
t83 = t126 + (-qJ(4) * t113 + t97) * t112;
t74 = t110 * t81 + t111 * t83;
t134 = t112 * t113;
t133 = t112 * t114;
t105 = -t114 * pkin(3) - pkin(2);
t92 = t110 * t112 - t111 * t114;
t93 = t110 * t114 + t111 * t112;
t78 = t92 * pkin(4) - t93 * qJ(5) + t105;
t131 = t78 * MDP(25);
t96 = pkin(3) * t134 + t113 * pkin(6);
t130 = t105 * MDP(21);
t129 = MDP(18) + MDP(22);
t128 = MDP(19) - MDP(24);
t125 = MDP(12) * t133;
t122 = t135 * t112;
t98 = t135 * t114;
t85 = t110 * t98 + t111 * t122;
t87 = -t110 * t122 + t111 * t98;
t90 = t93 * t113;
t91 = -t110 * t134 + t111 * t132;
t124 = t85 * t91 - t87 * t90;
t73 = -t110 * t83 + t111 * t81;
t103 = t111 * pkin(3) + pkin(4);
t121 = -t103 * MDP(25) - MDP(22);
t120 = t114 * MDP(13) - t112 * MDP(14);
t117 = t111 * MDP(18) - t110 * MDP(19);
t109 = t114 ^ 2;
t108 = t113 ^ 2;
t107 = t112 ^ 2;
t101 = t110 * pkin(3) + qJ(5);
t89 = t112 * t97 + t126;
t88 = -t112 * t136 + t94;
t75 = t90 * pkin(4) - t91 * qJ(5) + t96;
t72 = t115 * pkin(4) - t73;
t71 = -t115 * qJ(5) + t74;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t113 * MDP(10) + (t73 ^ 2 + t74 ^ 2 + t96 ^ 2) * MDP(21) + (t71 ^ 2 + t72 ^ 2 + t75 ^ 2) * MDP(25) + (t109 * MDP(11) + MDP(4) - 0.2e1 * t125) * t108 + (t115 * MDP(15) + 0.2e1 * pkin(1) * MDP(9) + 0.2e1 * (MDP(5) - t120) * t113) * t115 + 0.2e1 * (t108 * t138 - t88 * t115) * MDP(16) + 0.2e1 * (t108 * t137 + t89 * t115) * MDP(17) + (-t73 * t115 + t96 * t90) * t144 + (t74 * t115 + t96 * t91) * t143 + (-t73 * t91 - t74 * t90) * t142 + (t72 * t115 + t75 * t90) * t141 + (-t71 * t90 + t72 * t91) * t140 + (-t71 * t115 - t75 * t91) * t139; (t105 * t90 + t96 * t92) * MDP(18) + (t105 * t91 + t96 * t93) * MDP(19) + (-t73 * t93 - t74 * t92 + t124) * MDP(20) + (t96 * t105 - t73 * t85 + t74 * t87) * MDP(21) + (t75 * t92 + t78 * t90) * MDP(22) + (-t71 * t92 + t72 * t93 + t124) * MDP(23) + (-t75 * t93 - t78 * t91) * MDP(24) + (t71 * t87 + t72 * t85 + t75 * t78) * MDP(25) + (-pkin(6) * MDP(10) + t128 * t87 + t129 * t85 + MDP(7) - t145) * t115 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t133 + (-t107 + t109) * MDP(12) + (-pkin(2) * t112 - t137) * MDP(16) + (-pkin(2) * t114 + t138) * MDP(17)) * t113; MDP(8) + t107 * MDP(11) + 0.2e1 * t125 + (-0.2e1 * t93 * MDP(24) + t92 * t141 + t131) * t78 + (t93 * t143 + t92 * t144 + t130) * t105 + 0.2e1 * (t114 * MDP(16) - t112 * MDP(17)) * pkin(2) + t146 * (t85 ^ 2 + t87 ^ 2) + (t142 + t140) * (t85 * t93 - t87 * t92); t88 * MDP(16) - t89 * MDP(17) + t73 * MDP(18) - t74 * MDP(19) + t73 * MDP(22) + (-t101 * t90 - t103 * t91) * MDP(23) + t74 * MDP(24) + (t71 * t101 - t72 * t103) * MDP(25) + t120 * t113 + (-MDP(15) + (-pkin(4) - t103) * MDP(22) + (-qJ(5) - t101) * MDP(24)) * t115 + ((-t110 * t90 - t111 * t91) * MDP(20) + (t110 * t74 + t111 * t73) * MDP(21) - t117 * t115) * pkin(3); (-t101 * t92 - t103 * t93) * MDP(23) + (MDP(25) * t101 - t128) * t87 + (-MDP(18) + t121) * t85 + ((-t110 * t92 - t111 * t93) * MDP(20) + (t110 * t87 - t111 * t85) * MDP(21)) * pkin(3) + t145; MDP(15) + (t101 ^ 2 + t103 ^ 2) * MDP(25) + t103 * t141 + t101 * t139 + (0.2e1 * t117 + (t110 ^ 2 + t111 ^ 2) * MDP(21) * pkin(3)) * pkin(3); t96 * MDP(21) + t75 * MDP(25) + t128 * t91 + t129 * t90; t128 * t93 + t129 * t92 + t130 + t131; 0; t146; t115 * MDP(22) + t91 * MDP(23) + t72 * MDP(25); t93 * MDP(23) + t85 * MDP(25); t121; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
