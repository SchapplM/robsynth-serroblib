% Calculate joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:44
% EndTime: 2019-12-05 15:19:45
% DurationCPUTime: 0.24s
% Computational Cost: add. (149->71), mult. (403->125), div. (0->0), fcn. (442->12), ass. (0->46)
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t91 = -(MDP(18) * t68 + MDP(19) * t71) * pkin(9) + t68 * MDP(15) + t71 * MDP(16);
t90 = pkin(8) * t68;
t89 = pkin(8) * t71;
t72 = cos(qJ(4));
t88 = pkin(8) * t72;
t63 = sin(pkin(6));
t70 = sin(qJ(3));
t87 = t63 * t70;
t73 = cos(qJ(3));
t86 = t63 * t73;
t65 = cos(pkin(11));
t66 = cos(pkin(6));
t85 = t65 * t66;
t69 = sin(qJ(4));
t84 = t68 * t69;
t83 = t68 * t71;
t82 = t69 * t71;
t81 = t72 * MDP(17);
t80 = MDP(14) * t83;
t79 = MDP(15) * t71 - MDP(16) * t68;
t77 = t71 * MDP(18) - t68 * MDP(19);
t75 = t72 * MDP(11) - t69 * MDP(12) + MDP(4);
t74 = -MDP(11) - t77;
t67 = cos(pkin(5));
t64 = sin(pkin(5));
t62 = sin(pkin(11));
t61 = t71 ^ 2;
t60 = t69 ^ 2;
t59 = t68 ^ 2;
t57 = -t72 * pkin(4) - t69 * pkin(9) - pkin(3);
t56 = t69 * t66 + t72 * t87;
t55 = -t72 * t66 + t69 * t87;
t54 = -t64 * t65 * t63 + t67 * t66;
t53 = t68 * t57 + t71 * t88;
t52 = t71 * t57 - t68 * t88;
t51 = t71 * t56 - t68 * t86;
t50 = -t68 * t56 - t71 * t86;
t49 = t67 * t87 + (t62 * t73 + t70 * t85) * t64;
t48 = -t67 * t86 + (t62 * t70 - t73 * t85) * t64;
t47 = t49 * t72 + t54 * t69;
t46 = t49 * t69 - t54 * t72;
t45 = t47 * t71 + t48 * t68;
t44 = -t47 * t68 + t48 * t71;
t1 = [MDP(1) + (t67 ^ 2 + (t62 ^ 2 + t65 ^ 2) * t64 ^ 2) * MDP(2); t67 * MDP(2); MDP(2); -t49 * MDP(5) + (-t44 * t72 + t46 * t84) * MDP(18) + (t45 * t72 + t46 * t82) * MDP(19) - t75 * t48; (-t50 * t72 + t55 * t84) * MDP(18) + (t51 * t72 + t55 * t82) * MDP(19) + (-t70 * MDP(5) + t75 * t73) * t63; MDP(3) + (0.2e1 * pkin(3) * MDP(11) + t81) * t72 + (t61 * MDP(13) + MDP(6) - 0.2e1 * t80) * t60 + 0.2e1 * (-t52 * t72 + t60 * t90) * MDP(18) + 0.2e1 * (t53 * t72 + t60 * t89) * MDP(19) + 0.2e1 * (-pkin(3) * MDP(12) + (MDP(7) - t79) * t72) * t69; -t47 * MDP(12) + t74 * t46; -t56 * MDP(12) + t74 * t55; (-pkin(8) * MDP(12) + MDP(9) - t91) * t72 + (MDP(8) - pkin(8) * MDP(11) + MDP(13) * t83 + (-t59 + t61) * MDP(14) + (-pkin(4) * t68 - t89) * MDP(18) + (-pkin(4) * t71 + t90) * MDP(19)) * t69; t59 * MDP(13) + 0.2e1 * pkin(4) * t77 + MDP(10) + 0.2e1 * t80; t44 * MDP(18) - t45 * MDP(19); t50 * MDP(18) - t51 * MDP(19); t52 * MDP(18) - t53 * MDP(19) + t79 * t69 - t81; t91; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
