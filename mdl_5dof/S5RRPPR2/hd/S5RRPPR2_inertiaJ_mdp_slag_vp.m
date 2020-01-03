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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:37
% EndTime: 2020-01-03 11:57:38
% DurationCPUTime: 0.19s
% Computational Cost: add. (255->74), mult. (461->105), div. (0->0), fcn. (351->8), ass. (0->55)
t100 = sin(pkin(9));
t98 = t100 ^ 2;
t102 = cos(pkin(9));
t99 = t102 ^ 2;
t125 = t98 + t99;
t101 = sin(pkin(8));
t93 = t101 * pkin(2) + qJ(4);
t126 = t125 * t93;
t104 = sin(qJ(5));
t106 = cos(qJ(5));
t137 = t106 * MDP(17) - t104 * MDP(18);
t136 = -0.2e1 * t102;
t135 = 2 * MDP(10);
t134 = 0.2e1 * MDP(17);
t133 = 0.2e1 * MDP(18);
t105 = sin(qJ(2));
t132 = pkin(1) * t105;
t103 = cos(pkin(8));
t131 = t103 * pkin(2);
t123 = t106 * t98;
t121 = t102 * t106;
t112 = -t102 * pkin(4) - t100 * pkin(7) - pkin(3);
t107 = cos(qJ(2));
t95 = t107 * pkin(1) + pkin(2);
t76 = -t101 * t132 + t103 * t95;
t68 = t112 - t76;
t77 = t101 * t95 + t103 * t132;
t74 = qJ(4) + t77;
t64 = t104 * t68 + t74 * t121;
t130 = t64 * t102 + t74 * t123;
t78 = t112 - t131;
t67 = t104 * t78 + t93 * t121;
t129 = t67 * t102 + t93 * t123;
t128 = t125 * t74;
t75 = -pkin(3) - t76;
t94 = -pkin(3) - t131;
t127 = t75 + t94;
t124 = t104 * t98;
t97 = t100 * MDP(9);
t122 = t102 * t104;
t120 = t75 * MDP(11);
t119 = t94 * MDP(11);
t116 = t104 * t100 * MDP(15);
t89 = t106 * t100 * MDP(14);
t115 = t125 * MDP(11);
t114 = MDP(8) * t136 + 0.2e1 * t97;
t113 = t106 ^ 2 * t98 * MDP(12) - 0.2e1 * t104 * MDP(13) * t123 + t99 * MDP(16) + 0.2e1 * t102 * t116 + t89 * t136 + MDP(4);
t111 = (t107 * MDP(5) - t105 * MDP(6)) * pkin(1);
t110 = -t102 * MDP(16) - t116 + t89;
t109 = t97 + (-MDP(8) - t137) * t102;
t79 = t93 * t124;
t69 = t74 * t124;
t66 = t106 * t78 - t93 * t122;
t63 = t106 * t68 - t74 * t122;
t1 = [MDP(1) + (t76 ^ 2 + t77 ^ 2) * MDP(7) + t74 ^ 2 * t115 + (t114 + t120) * t75 + 0.2e1 * t111 + t128 * t135 + (-t63 * t102 + t69) * t134 + t130 * t133 + t113; (t126 + t128) * MDP(10) + (t74 * t126 + t75 * t94) * MDP(11) + (t69 + t79) * MDP(17) + (t129 + t130) * MDP(18) + t127 * t97 + (t101 * t77 + t103 * t76) * MDP(7) * pkin(2) + t111 + (-t127 * MDP(8) + (-t63 - t66) * MDP(17)) * t102 + t113; t93 ^ 2 * t115 + (t114 + t119) * t94 + (t101 ^ 2 + t103 ^ 2) * MDP(7) * pkin(2) ^ 2 + t126 * t135 + (-t66 * t102 + t79) * t134 + t129 * t133 + t113; 0; 0; MDP(7) + t115; t109 + t120; t109 + t119; 0; MDP(11); t63 * MDP(17) - t64 * MDP(18) + t110; t66 * MDP(17) - t67 * MDP(18) + t110; (-MDP(17) * t104 - MDP(18) * t106) * t100; t137; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
