% Calculate joint inertia matrix for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:51
% EndTime: 2019-12-05 16:15:52
% DurationCPUTime: 0.16s
% Computational Cost: add. (131->51), mult. (240->64), div. (0->0), fcn. (195->6), ass. (0->34)
t64 = sin(pkin(9));
t65 = cos(pkin(9));
t81 = t64 ^ 2 + t65 ^ 2;
t82 = t81 * qJ(4);
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t46 = t66 * t64 - t68 * t65;
t47 = t68 * t64 + t66 * t65;
t73 = t46 * MDP(17) + t47 * MDP(18);
t89 = t65 * MDP(8) - t64 * MDP(9);
t88 = 2 * MDP(10);
t69 = cos(qJ(3));
t87 = t69 * pkin(2);
t85 = t47 * MDP(14) - t46 * MDP(15);
t67 = sin(qJ(3));
t55 = t67 * pkin(2) + qJ(4);
t83 = t81 * t55;
t80 = pkin(3) * MDP(11);
t57 = -pkin(3) - t87;
t78 = t57 * MDP(11);
t77 = MDP(5) + (MDP(12) * t47 - 0.2e1 * MDP(13) * t46) * t47;
t56 = -t65 * pkin(4) - pkin(3);
t76 = t81 * MDP(11);
t75 = 0.2e1 * t89;
t74 = -t89 + t73;
t72 = (t69 * MDP(6) - t67 * MDP(7)) * pkin(2);
t71 = 0.2e1 * t73;
t61 = t65 * pkin(7);
t50 = t65 * qJ(4) + t61;
t49 = (-pkin(7) - qJ(4)) * t64;
t48 = t56 - t87;
t45 = t65 * t55 + t61;
t44 = (-pkin(7) - t55) * t64;
t1 = [MDP(1) + t76; 0; MDP(2) + t83 * t88 + t55 ^ 2 * t76 + (-t75 + t78) * t57 + t48 * t71 + 0.2e1 * t72 + t77; 0; (t82 + t83) * MDP(10) + (-t57 * pkin(3) + t55 * t82) * MDP(11) + t72 + t77 + t89 * (pkin(3) - t57) + t73 * (t48 + t56); t82 * t88 + qJ(4) ^ 2 * t76 + t56 * t71 + (t75 + t80) * pkin(3) + t77; 0; t74 + t78; t74 - t80; MDP(11); -t73; (t68 * t44 - t66 * t45) * MDP(17) + (-t66 * t44 - t68 * t45) * MDP(18) + t85; (t68 * t49 - t66 * t50) * MDP(17) + (-t66 * t49 - t68 * t50) * MDP(18) + t85; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
