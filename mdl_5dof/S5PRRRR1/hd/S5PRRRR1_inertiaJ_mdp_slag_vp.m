% Calculate joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaJ_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR1_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:24
% EndTime: 2019-12-05 17:03:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->50), mult. (333->77), div. (0->0), fcn. (349->8), ass. (0->33)
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t76 = -MDP(24) * t69 + MDP(25) * t65;
t74 = MDP(17) - t76;
t96 = pkin(2) * (MDP(24) * t65 + MDP(25) * t69);
t89 = t65 * MDP(21) + t69 * MDP(22);
t71 = cos(qJ(3));
t92 = pkin(2) * t71;
t95 = -0.2e1 * t92;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t93 = (-MDP(18) * t66 + t74 * t70) * pkin(2);
t67 = sin(qJ(3));
t57 = t66 * t71 + t70 * t67;
t68 = sin(qJ(2));
t49 = t57 * t68;
t91 = t49 * t57;
t90 = t65 * t69;
t56 = t66 * t67 - t70 * t71;
t83 = t56 * MDP(23);
t82 = t57 * MDP(18);
t81 = MDP(20) * t90;
t63 = t65 ^ 2;
t80 = t63 * MDP(19) + MDP(16) + 0.2e1 * t81;
t64 = t69 ^ 2;
t78 = (MDP(14) + MDP(19) * t90 + (-t63 + t64) * MDP(20)) * t57 + (-MDP(15) + t89) * t56;
t77 = MDP(21) * t69 - MDP(22) * t65;
t50 = t56 * t68;
t75 = t50 * MDP(18) - t74 * t49;
t72 = cos(qJ(2));
t44 = -t69 * t50 - t72 * t65;
t43 = t65 * t50 - t72 * t69;
t1 = [MDP(1); -t68 * MDP(4) + (t43 * t56 + t65 * t91) * MDP(24) + (-t44 * t56 + t69 * t91) * MDP(25) + (t71 * MDP(10) - t67 * MDP(11) - t56 * MDP(17) + MDP(3) - t82) * t72; t82 * t95 + MDP(2) + (MDP(5) * t67 + 0.2e1 * MDP(6) * t71) * t67 + (MDP(19) * t64 + MDP(12) - 0.2e1 * t81) * t57 ^ 2 + (t83 + 0.2e1 * (-MDP(13) + t77) * t57 + t74 * t95) * t56; (-MDP(10) * t67 - MDP(11) * t71) * t68 + t75; t67 * MDP(7) + t71 * MDP(8) + t78 + (-t56 * t66 - t57 * t70) * t96; MDP(9) + t80 + 0.2e1 * t93; t75; t78; t93 + t80; t80; t43 * MDP(24) - t44 * MDP(25); t77 * t57 + t76 * t92 + t83; -t66 * t96 + t89; t89; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
