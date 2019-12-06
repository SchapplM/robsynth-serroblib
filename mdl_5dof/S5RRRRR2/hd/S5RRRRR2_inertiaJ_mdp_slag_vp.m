% Calculate joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaJ_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:45
% EndTime: 2019-12-05 18:53:45
% DurationCPUTime: 0.23s
% Computational Cost: add. (255->71), mult. (664->102), div. (0->0), fcn. (649->8), ass. (0->46)
t108 = sin(qJ(4));
t109 = sin(qJ(3));
t112 = cos(qJ(4));
t113 = cos(qJ(3));
t93 = t108 * t109 - t112 * t113;
t145 = t93 * MDP(25);
t107 = sin(qJ(5));
t111 = cos(qJ(5));
t131 = t107 * MDP(23) + t111 * MDP(24);
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t144 = (-t110 * MDP(6) + (MDP(12) * t113 - MDP(13) * t109 + MDP(5)) * t114) * pkin(1);
t140 = pkin(1) * t110;
t94 = t108 * t113 + t112 * t109;
t80 = t94 * t140;
t81 = t93 * t140;
t142 = -t80 * MDP(19) + t81 * MDP(20);
t121 = -t111 * MDP(26) + t107 * MDP(27);
t119 = MDP(19) - t121;
t141 = (-MDP(20) * t108 + t119 * t112) * pkin(2);
t139 = t113 * pkin(2);
t138 = t108 * t93;
t137 = t112 * t94;
t136 = t80 * t107;
t135 = t80 * t111;
t134 = t107 * t111;
t132 = t94 * MDP(20);
t129 = t111 * t137;
t105 = t107 ^ 2;
t126 = MDP(22) * t134;
t128 = t105 * MDP(21) + MDP(18) + 0.2e1 * t126;
t127 = t107 * t94 * MDP(24);
t86 = t111 * t94 * MDP(23);
t95 = -t114 * pkin(1) - t139;
t68 = t107 * t81 + t111 * t95;
t125 = (t94 * t136 + t68 * t93) * MDP(26);
t69 = t107 * t95 - t111 * t81;
t124 = (t94 * t135 - t69 * t93) * MDP(27);
t106 = t111 ^ 2;
t123 = (MDP(16) + MDP(21) * t134 + (-t105 + t106) * MDP(22)) * t94 + (-MDP(17) + t131) * t93;
t122 = -t127 + t86 + t145;
t120 = (-t137 - t138) * t107;
t118 = t113 * MDP(10) + t109 * MDP(9) + t123;
t115 = MDP(4) + (0.2e1 * t86 + t145) * t93 - 0.2e1 * (MDP(15) * t94 + t127) * t93 + (MDP(7) * t109 + 0.2e1 * MDP(8) * t113) * t109 + (MDP(21) * t106 + MDP(14) - 0.2e1 * t126) * t94 ^ 2;
t82 = t111 * pkin(2) * t138;
t1 = [MDP(1) + t115 + 0.2e1 * t124 + 0.2e1 * t125 + 0.2e1 * (t93 * MDP(19) + t132) * t95 + 0.2e1 * t144; t115 + t144 + (t95 * MDP(19) - t119 * t139) * t93 + t124 + t125 + (t95 - t139) * t132; 0.2e1 * (-t119 * t93 - t132) * t139 + t115; (pkin(2) * t120 - t135) * MDP(26) + (-pkin(2) * t129 + t136 - t82) * MDP(27) + (-MDP(12) * t109 - MDP(13) * t113) * t140 + t118 + t142; -t82 * MDP(27) + (MDP(26) * t120 - MDP(27) * t129) * pkin(2) + t118; MDP(11) + t128 + 0.2e1 * t141; t121 * t80 + t123 + t142; t123; t141 + t128; t128; t68 * MDP(26) - t69 * MDP(27) + t122; t121 * t139 + t122; (-MDP(26) * t107 - MDP(27) * t111) * t108 * pkin(2) + t131; t131; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
