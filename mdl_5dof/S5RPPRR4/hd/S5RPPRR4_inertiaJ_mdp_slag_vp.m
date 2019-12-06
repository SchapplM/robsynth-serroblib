% Calculate joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR4_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:56
% EndTime: 2019-12-05 17:44:58
% DurationCPUTime: 0.30s
% Computational Cost: add. (378->78), mult. (795->115), div. (0->0), fcn. (811->8), ass. (0->43)
t101 = cos(qJ(5));
t100 = sin(qJ(4));
t102 = cos(qJ(4));
t96 = sin(pkin(8));
t98 = cos(pkin(8));
t86 = -t98 * pkin(2) - t96 * qJ(3) - pkin(1);
t97 = cos(pkin(9));
t82 = t97 * t86;
t95 = sin(pkin(9));
t71 = -t97 * t96 * pkin(6) + t82 + (-qJ(2) * t95 - pkin(3)) * t98;
t122 = t95 * t96;
t117 = qJ(2) * t98;
t77 = t97 * t117 + t95 * t86;
t75 = -pkin(6) * t122 + t77;
t62 = -t100 * t75 + t102 * t71;
t83 = -t100 * t95 + t102 * t97;
t80 = t83 * t96;
t60 = -t98 * pkin(4) - t80 * pkin(7) + t62;
t63 = t100 * t71 + t102 * t75;
t84 = t100 * t97 + t102 * t95;
t79 = t84 * t96;
t61 = -t79 * pkin(7) + t63;
t99 = sin(qJ(5));
t66 = t101 * t79 + t99 * t80;
t67 = t101 * t80 - t99 * t79;
t127 = -(t101 * t61 + t99 * t60) * MDP(25) + (t101 * t60 - t99 * t61) * MDP(24) + t67 * MDP(21) - t66 * MDP(22);
t126 = t80 * MDP(14) - t79 * MDP(15) + t62 * MDP(17) - t63 * MDP(18) + t127;
t104 = (MDP(24) * t101 - MDP(25) * t99) * pkin(4);
t120 = (t101 * t83 - t99 * t84) * MDP(24) - (t101 * t84 + t99 * t83) * MDP(25);
t124 = t83 * MDP(17) - t84 * MDP(18) + t120;
t119 = pkin(3) * t122 + t96 * qJ(2);
t118 = t95 ^ 2 + t97 ^ 2;
t114 = MDP(16) + MDP(23);
t76 = -t95 * t117 + t82;
t113 = t76 * t97 + t77 * t95;
t112 = t95 * MDP(8) + t97 * MDP(9);
t109 = t79 * MDP(17) + t80 * MDP(18);
t106 = t66 * MDP(24) + t67 * MDP(25);
t103 = qJ(2) ^ 2;
t94 = t98 ^ 2;
t92 = t96 ^ 2;
t89 = t92 * t103;
t1 = [MDP(1) + (pkin(1) ^ 2 + t89) * MDP(7) + (t76 ^ 2 + t77 ^ 2 + t89) * MDP(11) + (t103 * MDP(7) + t114) * t94 + (MDP(12) * t80 - 0.2e1 * t79 * MDP(13)) * t80 + (MDP(19) * t67 - 0.2e1 * t66 * MDP(20)) * t67 + 0.2e1 * t109 * t119 + 0.2e1 * t106 * (t79 * pkin(4) + t119) + 0.2e1 * (t94 * MDP(6) + (MDP(6) + t112) * t92) * qJ(2) + 0.2e1 * (-t113 * MDP(10) - pkin(1) * MDP(5)) * t96 + 0.2e1 * (pkin(1) * MDP(4) - t76 * MDP(8) + t77 * MDP(9) - t126) * t98; -pkin(1) * MDP(7) + t113 * MDP(11) + (-t118 * MDP(10) + MDP(5)) * t96 + (-t97 * MDP(8) + t95 * MDP(9) - MDP(4) - t124) * t98; t118 * MDP(11) + MDP(7); (MDP(11) * qJ(2) + t112) * t96 + t106 + t109; 0; MDP(11); (-t114 - t104) * t98 + t126; t124; 0; 0.2e1 * t104 + t114; -t98 * MDP(23) + t127; t120; 0; MDP(23) + t104; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
