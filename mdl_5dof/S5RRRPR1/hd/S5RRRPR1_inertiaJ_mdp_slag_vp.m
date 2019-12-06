% Calculate joint inertia matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:52
% EndTime: 2019-12-05 18:38:52
% DurationCPUTime: 0.22s
% Computational Cost: add. (564->92), mult. (1034->139), div. (0->0), fcn. (1157->8), ass. (0->51)
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t103 = sin(qJ(3));
t104 = sin(qJ(2));
t106 = cos(qJ(3));
t107 = cos(qJ(2));
t89 = t103 * t104 - t106 * t107;
t90 = t103 * t107 + t104 * t106;
t76 = -t100 * t90 - t101 * t89;
t99 = -pkin(2) * t107 - pkin(1);
t81 = pkin(3) * t89 + t99;
t126 = -0.2e1 * pkin(4) * t76 + 0.2e1 * t81;
t125 = 0.2e1 * t99;
t124 = 0.2e1 * t107;
t123 = pkin(6) + pkin(7);
t122 = pkin(2) * t103;
t121 = pkin(3) * t100;
t94 = t123 * t104;
t95 = t123 * t107;
t113 = -t103 * t95 - t106 * t94;
t73 = -qJ(4) * t90 + t113;
t111 = t103 * t94 - t106 * t95;
t74 = -qJ(4) * t89 - t111;
t62 = t100 * t73 + t101 * t74;
t102 = sin(qJ(5));
t96 = pkin(3) * t101 + pkin(4);
t120 = t102 * t96;
t105 = cos(qJ(5));
t77 = -t100 * t89 + t101 * t90;
t65 = t102 * t77 - t105 * t76;
t119 = t65 * MDP(25);
t98 = pkin(2) * t106 + pkin(3);
t83 = -t100 * t122 + t101 * t98;
t82 = pkin(4) + t83;
t85 = t100 * t98 + t101 * t122;
t71 = -t102 * t85 + t105 * t82;
t118 = t71 * MDP(25);
t72 = -t102 * t82 - t105 * t85;
t117 = t72 * MDP(26);
t93 = t105 * t96;
t116 = (-t102 * t121 + t93) * MDP(25);
t115 = (-t105 * t121 - t120) * MDP(26);
t114 = MDP(15) + MDP(24);
t61 = -t100 * t74 + t101 * t73;
t59 = -pkin(8) * t77 + t61;
t60 = pkin(8) * t76 + t62;
t66 = t102 * t76 + t105 * t77;
t112 = t66 * MDP(22) - t65 * MDP(23) + (-t102 * t60 + t105 * t59) * MDP(25) + (-t102 * t59 - t105 * t60) * MDP(26);
t110 = (MDP(16) * t106 - MDP(17) * t103) * pkin(2);
t109 = t90 * MDP(13) - t89 * MDP(14) + t113 * MDP(16) + t111 * MDP(17) + t112;
t1 = [MDP(1) + pkin(1) * MDP(9) * t124 + t89 * MDP(16) * t125 + 0.2e1 * (-t61 * t77 + t62 * t76) * MDP(18) + (t61 ^ 2 + t62 ^ 2 + t81 ^ 2) * MDP(19) + t119 * t126 + (MDP(11) * t90 - 0.2e1 * MDP(12) * t89 + MDP(17) * t125) * t90 + (MDP(20) * t66 - 0.2e1 * MDP(21) * t65 + MDP(26) * t126) * t66 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t104 + MDP(5) * t124) * t104; t104 * MDP(6) + t107 * MDP(7) + (t76 * t85 - t77 * t83) * MDP(18) + (t61 * t83 + t62 * t85) * MDP(19) + (-MDP(10) * t107 - MDP(9) * t104) * pkin(6) + t109; MDP(8) + (t83 ^ 2 + t85 ^ 2) * MDP(19) + 0.2e1 * t110 + 0.2e1 * t118 + 0.2e1 * t117 + t114; ((t100 * t76 - t101 * t77) * MDP(18) + (t100 * t62 + t101 * t61) * MDP(19)) * pkin(3) + t109; (t71 + t93) * MDP(25) + (t72 - t120) * MDP(26) + (t101 * t83 * MDP(19) + (MDP(19) * t85 - MDP(25) * t102 - MDP(26) * t105) * t100) * pkin(3) + t110 + t114; (t100 ^ 2 + t101 ^ 2) * MDP(19) * pkin(3) ^ 2 + 0.2e1 * t116 + 0.2e1 * t115 + t114; MDP(19) * t81 + MDP(26) * t66 + t119; 0; 0; MDP(19); t112; MDP(24) + t117 + t118; MDP(24) + t115 + t116; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
