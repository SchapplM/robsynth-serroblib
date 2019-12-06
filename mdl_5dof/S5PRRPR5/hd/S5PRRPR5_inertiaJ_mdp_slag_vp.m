% Calculate joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:56
% EndTime: 2019-12-05 16:27:58
% DurationCPUTime: 0.23s
% Computational Cost: add. (246->83), mult. (535->144), div. (0->0), fcn. (580->10), ass. (0->44)
t72 = sin(pkin(10));
t74 = cos(pkin(10));
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t63 = t72 * t79 + t74 * t76;
t75 = sin(qJ(5));
t98 = t63 * t75;
t78 = cos(qJ(5));
t97 = t63 * t78;
t73 = sin(pkin(5));
t77 = sin(qJ(2));
t96 = t73 * t77;
t80 = cos(qJ(2));
t95 = t73 * t80;
t94 = t75 * t78;
t93 = -qJ(4) - pkin(7);
t92 = cos(pkin(5));
t91 = MDP(10) * t79;
t62 = t72 * t76 - t74 * t79;
t90 = t62 * MDP(18);
t69 = -t79 * pkin(3) - pkin(2);
t89 = t69 * MDP(13);
t88 = MDP(15) * t94;
t87 = t93 * t76;
t86 = t78 * MDP(19) - t75 * MDP(20);
t85 = MDP(19) * t75 + MDP(20) * t78;
t84 = (MDP(16) * t78 - MDP(17) * t75) * t63;
t83 = -t76 * t96 + t92 * t79;
t82 = t75 * MDP(16) + t78 * MDP(17) - t85 * (t72 * pkin(3) + pkin(8));
t71 = t78 ^ 2;
t70 = t75 ^ 2;
t68 = -t74 * pkin(3) - pkin(4);
t65 = t93 * t79;
t60 = t92 * t76 + t79 * t96;
t58 = -t74 * t65 + t72 * t87;
t56 = -t72 * t65 - t74 * t87;
t55 = t62 * pkin(4) - t63 * pkin(8) + t69;
t54 = t74 * t60 + t72 * t83;
t52 = t72 * t60 - t74 * t83;
t51 = t78 * t54 - t75 * t95;
t50 = -t75 * t54 - t78 * t95;
t49 = t75 * t55 + t78 * t58;
t48 = t78 * t55 - t75 * t58;
t1 = [MDP(1) + (t73 ^ 2 * t80 ^ 2 + t52 ^ 2 + t54 ^ 2) * MDP(13); (t52 * t63 - t54 * t62) * MDP(12) + (t52 * t56 + t54 * t58) * MDP(13) + (t50 * t62 + t52 * t98) * MDP(19) + (-t51 * t62 + t52 * t97) * MDP(20) + (-t77 * MDP(4) + (-MDP(11) * t76 + MDP(3) - t89 + t91) * t80) * t73; MDP(2) + 0.2e1 * pkin(2) * t91 + (t56 ^ 2 + t58 ^ 2 + t69 ^ 2) * MDP(13) + (t71 * MDP(14) - 0.2e1 * t88) * t63 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t76 + 0.2e1 * t79 * MDP(6)) * t76 + (0.2e1 * t84 + t90) * t62 + 0.2e1 * (t56 * t63 - t58 * t62) * MDP(12) + 0.2e1 * (t48 * t62 + t56 * t98) * MDP(19) + 0.2e1 * (-t49 * t62 + t56 * t97) * MDP(20); t83 * MDP(10) - t60 * MDP(11) - t86 * t52 + (-t52 * t74 + t54 * t72) * MDP(13) * pkin(3); t76 * MDP(7) + t79 * MDP(8) - t86 * t56 + t82 * t62 + (-t76 * MDP(10) - t79 * MDP(11)) * pkin(7) + (MDP(14) * t94 + (-t70 + t71) * MDP(15) + t85 * t68) * t63 + ((-t62 * t72 - t63 * t74) * MDP(12) + (-t56 * t74 + t58 * t72) * MDP(13)) * pkin(3); 0.2e1 * t88 + t70 * MDP(14) + MDP(9) + (t72 ^ 2 + t74 ^ 2) * MDP(13) * pkin(3) ^ 2 - 0.2e1 * t86 * t68; -MDP(13) * t95; t62 * t86 + t89; 0; MDP(13); t50 * MDP(19) - t51 * MDP(20); t48 * MDP(19) - t49 * MDP(20) + t84 + t90; t82; t86; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
