% Calculate joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:52
% EndTime: 2022-01-20 10:05:53
% DurationCPUTime: 0.19s
% Computational Cost: add. (241->71), mult. (440->106), div. (0->0), fcn. (335->8), ass. (0->50)
t97 = sin(pkin(9));
t95 = t97 ^ 2;
t99 = cos(pkin(9));
t96 = t99 ^ 2;
t119 = t95 + t96;
t98 = sin(pkin(8));
t91 = t98 * pkin(2) + qJ(4);
t120 = t119 * t91;
t101 = sin(qJ(5));
t103 = cos(qJ(5));
t130 = t103 * MDP(16) - t101 * MDP(17);
t129 = -0.2e1 * t99;
t128 = 2 * MDP(9);
t127 = 0.2e1 * MDP(16);
t126 = 0.2e1 * MDP(17);
t102 = sin(qJ(2));
t125 = pkin(1) * t102;
t100 = cos(pkin(8));
t124 = t100 * pkin(2);
t116 = t103 * t95;
t115 = t103 * t99;
t110 = -t99 * pkin(4) - t97 * pkin(7) - pkin(3);
t104 = cos(qJ(2));
t93 = t104 * pkin(1) + pkin(2);
t74 = t100 * t93 - t98 * t125;
t66 = t110 - t74;
t75 = t100 * t125 + t98 * t93;
t72 = qJ(4) + t75;
t62 = t101 * t66 + t72 * t115;
t123 = t72 * t116 + t62 * t99;
t76 = t110 - t124;
t65 = t101 * t76 + t91 * t115;
t122 = t91 * t116 + t65 * t99;
t121 = t119 * t72;
t118 = t101 * t95;
t117 = t101 * t99;
t112 = MDP(8) * t129;
t111 = t101 * t97 * MDP(14);
t87 = t103 * t97 * MDP(13);
t109 = t103 ^ 2 * t95 * MDP(11) - 0.2e1 * t101 * MDP(12) * t116 + t96 * MDP(15) + 0.2e1 * t99 * t111 + t87 * t129 + MDP(4);
t108 = (t104 * MDP(5) - t102 * MDP(6)) * pkin(1);
t107 = (-MDP(8) - t130) * t99;
t106 = -t99 * MDP(15) - t111 + t87;
t92 = -pkin(3) - t124;
t77 = t91 * t118;
t73 = -pkin(3) - t74;
t67 = t72 * t118;
t64 = t103 * t76 - t91 * t117;
t61 = t103 * t66 - t72 * t117;
t1 = [MDP(1) + (t74 ^ 2 + t75 ^ 2) * MDP(7) + t73 * t112 + (t119 * t72 ^ 2 + t73 ^ 2) * MDP(10) + 0.2e1 * t108 + t121 * t128 + (-t61 * t99 + t67) * t127 + t123 * t126 + t109; (t120 + t121) * MDP(9) + (t72 * t120 + t73 * t92) * MDP(10) + (t67 + t77) * MDP(16) + (t122 + t123) * MDP(17) + (t100 * t74 + t75 * t98) * MDP(7) * pkin(2) + t108 + ((-t73 - t92) * MDP(8) + (-t61 - t64) * MDP(16)) * t99 + t109; t92 * t112 + (t119 * t91 ^ 2 + t92 ^ 2) * MDP(10) + (t100 ^ 2 + t98 ^ 2) * MDP(7) * pkin(2) ^ 2 + t120 * t128 + (-t64 * t99 + t77) * t127 + t122 * t126 + t109; 0; 0; t119 * MDP(10) + MDP(7); t73 * MDP(10) + t107; t92 * MDP(10) + t107; 0; MDP(10); t61 * MDP(16) - t62 * MDP(17) + t106; t64 * MDP(16) - t65 * MDP(17) + t106; (-MDP(16) * t101 - MDP(17) * t103) * t97; t130; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
