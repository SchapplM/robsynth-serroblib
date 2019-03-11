% Calculate joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPPRR2_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:02
% EndTime: 2019-03-09 01:32:03
% DurationCPUTime: 0.22s
% Computational Cost: add. (261->70), mult. (445->102), div. (0->0), fcn. (441->8), ass. (0->43)
t83 = sin(qJ(6));
t85 = cos(qJ(6));
t89 = t83 * MDP(24) + t85 * MDP(25);
t79 = sin(pkin(10));
t81 = cos(pkin(10));
t99 = t79 ^ 2 + t81 ^ 2;
t108 = MDP(11) * t99;
t82 = cos(pkin(9));
t72 = -t82 * pkin(1) - pkin(2);
t107 = t72 * MDP(7);
t102 = cos(qJ(5));
t84 = sin(qJ(5));
t61 = -t102 * t81 + t84 * t79;
t62 = t102 * t79 + t84 * t81;
t90 = t85 * MDP(24) - t83 * MDP(25);
t88 = MDP(17) + t90;
t105 = t61 * MDP(18) - t62 * t88;
t104 = t61 ^ 2;
t59 = t62 ^ 2;
t68 = -qJ(4) + t72;
t103 = -pkin(7) + t68;
t57 = t103 * t79;
t58 = t103 * t81;
t54 = -t102 * t58 + t84 * t57;
t101 = t54 * t61;
t100 = t83 * t85;
t98 = t79 * MDP(8);
t97 = t81 * MDP(9);
t96 = MDP(7) + t108;
t80 = sin(pkin(9));
t69 = t80 * pkin(1) + qJ(3);
t64 = t79 * pkin(4) + t69;
t92 = MDP(20) * t100;
t91 = -MDP(21) * t85 + MDP(22) * t83;
t87 = t83 * MDP(21) + t85 * MDP(22) - t89 * pkin(8);
t78 = t85 ^ 2;
t77 = t83 ^ 2;
t56 = t99 * t68;
t55 = t102 * t57 + t84 * t58;
t53 = t62 * pkin(5) + t61 * pkin(8) + t64;
t52 = t83 * t53 + t85 * t55;
t51 = t85 * t53 - t83 * t55;
t1 = [0.2e1 * t64 * t62 * MDP(17) + t59 * MDP(23) + MDP(1) + (t80 ^ 2 + t82 ^ 2) * MDP(4) * pkin(1) ^ 2 + t68 ^ 2 * t108 + (t78 * MDP(19) + MDP(12) - 0.2e1 * t92) * t104 + 0.2e1 * (-t64 * MDP(18) + (MDP(13) + t91) * t62) * t61 - 0.2e1 * t56 * MDP(10) + 0.2e1 * (-t83 * t101 + t51 * t62) * MDP(24) + 0.2e1 * (-t85 * t101 - t52 * t62) * MDP(25) + ((2 * MDP(5)) + t107) * t72 + (0.2e1 * t98 + 0.2e1 * t97 + (2 * MDP(6)) + (MDP(7) + MDP(11)) * t69) * t69; 0; MDP(4) + t96; -t99 * MDP(10) + t56 * MDP(11) + MDP(5) + t107 + t89 * (-t59 - t104); 0; t96; t69 * MDP(11) - t105 + t97 + t98; 0; 0; MDP(11); -t55 * MDP(18) - t88 * t54 + (-MDP(15) + t87) * t62 + (-MDP(14) - MDP(19) * t100 + (t77 - t78) * MDP(20) + t89 * pkin(5)) * t61; t105; -t62 * MDP(18) - t61 * t88; 0; t77 * MDP(19) + 0.2e1 * pkin(5) * t90 + MDP(16) + 0.2e1 * t92; t62 * MDP(23) + t51 * MDP(24) - t52 * MDP(25) + t91 * t61; t89 * t61; -t89 * t62; t90; t87; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
