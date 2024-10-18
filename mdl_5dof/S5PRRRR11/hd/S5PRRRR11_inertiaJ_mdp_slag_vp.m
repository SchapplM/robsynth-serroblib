% Calculate joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR11_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:58
% DurationCPUTime: 0.27s
% Computational Cost: add. (380->96), mult. (955->135), div. (0->0), fcn. (984->8), ass. (0->58)
t91 = sin(pkin(5));
t98 = cos(qJ(3));
t114 = t91 * t98;
t95 = sin(qJ(3));
t115 = t91 * t95;
t94 = sin(qJ(4));
t97 = cos(qJ(4));
t76 = -t97 * t114 + t94 * t115;
t77 = (t94 * t98 + t95 * t97) * t91;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t65 = t96 * t76 + t93 * t77;
t66 = -t93 * t76 + t96 * t77;
t125 = t66 * MDP(21) - t65 * MDP(22);
t124 = t77 * MDP(14) - t76 * MDP(15);
t123 = pkin(7) + pkin(8);
t122 = pkin(2) * t98;
t121 = pkin(3) * t94;
t92 = cos(pkin(5));
t120 = t92 * pkin(3);
t119 = t92 * pkin(4);
t118 = t96 * pkin(4);
t117 = t97 * pkin(3);
t89 = t91 ^ 2;
t116 = t89 * t95;
t105 = t92 * t95 * pkin(2);
t70 = t123 * t114 + t105;
t112 = t97 * t70;
t85 = t92 * t122;
t69 = -t123 * t115 + t120 + t85;
t60 = t94 * t69 + t112;
t58 = -t76 * pkin(9) + t60;
t113 = t96 * t58;
t111 = -t65 * MDP(24) - t66 * MDP(25);
t86 = pkin(4) + t117;
t83 = t96 * t86;
t78 = -t93 * t121 + t83;
t110 = t78 * MDP(24);
t79 = t96 * t121 + t93 * t86;
t109 = t79 * MDP(25);
t108 = t93 * MDP(25);
t107 = t97 * MDP(17);
t106 = MDP(16) + MDP(23);
t104 = t92 * MDP(23) + t125;
t59 = t97 * t69 - t94 * t70;
t57 = -t77 * pkin(9) + t119 + t59;
t54 = t96 * t57 - t93 * t58;
t103 = MDP(9) + t106;
t102 = -t76 * MDP(17) - t77 * MDP(18) + t111;
t82 = (-pkin(3) * t98 - pkin(2)) * t91;
t55 = t93 * t57 + t113;
t101 = t92 * MDP(16) + t104 + t124;
t100 = (t95 * MDP(7) + t98 * MDP(8)) * t91;
t99 = (t96 * MDP(24) - t108) * pkin(4);
t81 = pkin(7) * t114 + t105;
t80 = -pkin(7) * t115 + t85;
t67 = t76 * pkin(4) + t82;
t1 = [MDP(1); 0; MDP(2) + (MDP(5) * t95 + 0.2e1 * MDP(6) * t98) * t116 + (MDP(12) * t77 - 0.2e1 * t76 * MDP(13)) * t77 + (MDP(19) * t66 - 0.2e1 * t65 * MDP(20)) * t66 + t103 * t92 ^ 2 + 0.2e1 * (t100 + t124 + t125) * t92 + 0.2e1 * (t89 * t122 + t80 * t92) * MDP(10) + 0.2e1 * (-pkin(2) * t116 - t81 * t92) * MDP(11) + 0.2e1 * (t59 * t92 + t82 * t76) * MDP(17) + 0.2e1 * (-t60 * t92 + t82 * t77) * MDP(18) + 0.2e1 * (t54 * t92 + t67 * t65) * MDP(24) + 0.2e1 * (-t55 * t92 + t67 * t66) * MDP(25); (MDP(10) * t98 - MDP(11) * t95) * t91 + t102; t92 * MDP(9) + t80 * MDP(10) - t81 * MDP(11) + (t92 * t117 + t59) * MDP(17) + (-t112 + (-t69 - t120) * t94) * MDP(18) + (t78 * t92 + t54) * MDP(24) + (-t79 * t92 - t55) * MDP(25) + t100 + t101; 0.2e1 * (-t94 * MDP(18) + t107) * pkin(3) + 0.2e1 * t110 - 0.2e1 * t109 + t103; t102; t59 * MDP(17) - t60 * MDP(18) + (t92 * t118 + t54) * MDP(24) + (-t113 + (-t57 - t119) * t93) * MDP(25) + t101; (t83 + t118) * MDP(24) + (-pkin(4) - t86) * t108 + (t107 + (-MDP(24) * t93 - MDP(25) * t96 - MDP(18)) * t94) * pkin(3) + t106; t106 + 0.2e1 * t99; t111; t54 * MDP(24) - t55 * MDP(25) + t104; MDP(23) - t109 + t110; MDP(23) + t99; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
