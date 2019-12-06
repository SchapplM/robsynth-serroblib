% Calculate joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR7_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:13:00
% EndTime: 2019-12-05 17:13:01
% DurationCPUTime: 0.19s
% Computational Cost: add. (239->71), mult. (473->102), div. (0->0), fcn. (507->8), ass. (0->41)
t78 = sin(qJ(4));
t79 = sin(qJ(3));
t82 = cos(qJ(4));
t83 = cos(qJ(3));
t68 = t78 * t79 - t82 * t83;
t76 = -pkin(3) * t83 - pkin(2);
t104 = 0.2e1 * pkin(4) * t68 + 0.2e1 * t76;
t103 = 0.2e1 * t76;
t102 = pkin(6) + pkin(7);
t101 = pkin(3) * t78;
t69 = t78 * t83 + t79 * t82;
t80 = sin(qJ(2));
t61 = t69 * t80;
t62 = t68 * t80;
t77 = sin(qJ(5));
t81 = cos(qJ(5));
t100 = (-t61 * t81 + t62 * t77) * MDP(24) + (t61 * t77 + t62 * t81) * MDP(25);
t99 = MDP(10) * t83;
t98 = MDP(17) * t68;
t54 = t68 * t81 + t69 * t77;
t97 = MDP(24) * t54;
t75 = pkin(3) * t82 + pkin(4);
t73 = t81 * t75;
t96 = (-t101 * t77 + t73) * MDP(24);
t95 = (-t101 * t81 - t75 * t77) * MDP(25);
t94 = t77 * MDP(25);
t93 = t82 * MDP(17);
t92 = MDP(16) + MDP(23);
t71 = t102 * t79;
t72 = t102 * t83;
t91 = -t82 * t71 - t72 * t78;
t50 = -pkin(8) * t69 + t91;
t88 = t71 * t78 - t72 * t82;
t51 = -pkin(8) * t68 - t88;
t55 = -t68 * t77 + t69 * t81;
t90 = t55 * MDP(21) - t54 * MDP(22) + (t50 * t81 - t51 * t77) * MDP(24) + (-t50 * t77 - t51 * t81) * MDP(25);
t89 = -t61 * MDP(17) + t62 * MDP(18) + t100;
t87 = -MDP(10) * t79 - MDP(11) * t83;
t86 = (MDP(24) * t81 - t94) * pkin(4);
t85 = t69 * MDP(14) - t68 * MDP(15) + t91 * MDP(17) + t88 * MDP(18) + t90;
t1 = [MDP(1); -t80 * MDP(4) + (-MDP(11) * t79 - MDP(18) * t69 - MDP(25) * t55 + MDP(3) - t97 - t98 + t99) * cos(qJ(2)); 0.2e1 * pkin(2) * t99 + t98 * t103 + t97 * t104 + MDP(2) + (-0.2e1 * MDP(11) * pkin(2) + MDP(5) * t79 + 0.2e1 * MDP(6) * t83) * t79 + (MDP(12) * t69 - 0.2e1 * MDP(13) * t68 + MDP(18) * t103) * t69 + (MDP(19) * t55 - 0.2e1 * MDP(20) * t54 + MDP(25) * t104) * t55; t80 * t87 + t89; t79 * MDP(7) + t83 * MDP(8) + pkin(6) * t87 + t85; MDP(9) + 0.2e1 * (-MDP(18) * t78 + t93) * pkin(3) + 0.2e1 * t96 + 0.2e1 * t95 + t92; t89; t85; (pkin(4) * t81 + t73) * MDP(24) + (-pkin(4) - t75) * t94 + (t93 + (-MDP(24) * t77 - MDP(25) * t81 - MDP(18)) * t78) * pkin(3) + t92; 0.2e1 * t86 + t92; t100; t90; MDP(23) + t95 + t96; MDP(23) + t86; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
