% Calculate joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP7_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:27
% DurationCPUTime: 0.28s
% Computational Cost: add. (206->98), mult. (454->148), div. (0->0), fcn. (429->8), ass. (0->45)
t85 = 2 * MDP(19);
t61 = sin(qJ(4));
t84 = pkin(7) * t61;
t64 = cos(qJ(4));
t83 = pkin(7) * t64;
t65 = cos(qJ(3));
t82 = pkin(7) * t65;
t59 = sin(pkin(5));
t63 = sin(qJ(2));
t81 = t59 * t63;
t66 = cos(qJ(2));
t80 = t59 * t66;
t79 = t61 * t64;
t78 = -qJ(5) - pkin(8);
t77 = MDP(20) * pkin(4);
t62 = sin(qJ(3));
t76 = qJ(5) * t62;
t60 = cos(pkin(5));
t49 = t60 * t62 + t65 * t81;
t44 = t49 * t64 - t61 * t80;
t75 = t44 * MDP(18);
t55 = -t64 * pkin(4) - pkin(3);
t74 = t55 * MDP(20);
t73 = t61 * MDP(15);
t72 = t65 * MDP(16);
t71 = t64 * t82;
t70 = MDP(13) * t79;
t69 = -MDP(19) * pkin(4) + MDP(14);
t68 = t64 * MDP(17) - t61 * MDP(18);
t67 = t61 * MDP(17) + t64 * MDP(18);
t58 = t64 ^ 2;
t57 = t62 ^ 2;
t56 = t61 ^ 2;
t54 = t78 * t64;
t53 = t78 * t61;
t52 = -t65 * pkin(3) - t62 * pkin(8) - pkin(2);
t51 = (pkin(4) * t61 + pkin(7)) * t62;
t50 = t64 * t52;
t48 = -t60 * t65 + t62 * t81;
t47 = t61 * t52 + t71;
t46 = -t61 * t82 + t50;
t45 = t71 + (t52 - t76) * t61;
t43 = -t49 * t61 - t64 * t80;
t42 = -t64 * t76 + t50 + (-pkin(4) - t84) * t65;
t1 = [MDP(1) + (t43 ^ 2 + t44 ^ 2 + t48 ^ 2) * MDP(20); (t43 * t42 + t44 * t45 + t48 * t51) * MDP(20) + (-t43 * MDP(17) + t75) * t65 + ((-t43 * t64 - t44 * t61) * MDP(19) + t67 * t48) * t62 + (-t63 * MDP(4) + (MDP(10) * t65 - MDP(11) * t62 + MDP(3)) * t66) * t59; MDP(2) + (t42 ^ 2 + t45 ^ 2 + t51 ^ 2) * MDP(20) + (0.2e1 * pkin(2) * MDP(10) + t72) * t65 + (t58 * MDP(12) + MDP(5) - 0.2e1 * t70) * t57 + 0.2e1 * (-t46 * t65 + t57 * t84) * MDP(17) + 0.2e1 * (t47 * t65 + t57 * t83) * MDP(18) + (-0.2e1 * pkin(2) * MDP(11) + (-t42 * t64 - t45 * t61) * t85 + 0.2e1 * (-t64 * MDP(14) + MDP(6) + t73) * t65) * t62; -t49 * MDP(11) + (-t43 * t61 + t44 * t64) * MDP(19) + (t43 * t53 - t44 * t54) * MDP(20) + (-MDP(10) - t68 + t74) * t48; (-t42 * t61 + t45 * t64) * MDP(19) + (t42 * t53 - t45 * t54 + t51 * t55) * MDP(20) + (-pkin(7) * MDP(11) - t61 * MDP(14) - t64 * MDP(15) + pkin(8) * t67 + MDP(8)) * t65 + (MDP(7) - pkin(7) * MDP(10) + MDP(12) * t79 + (-t56 + t58) * MDP(13) + (-pkin(3) * t61 - t83) * MDP(17) + (-pkin(3) * t64 + t84) * MDP(18) + (-t53 * t64 + t54 * t61) * MDP(19)) * t62; MDP(9) + t56 * MDP(12) + 0.2e1 * t70 + (-t53 * t61 - t54 * t64) * t85 + (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) * MDP(20) + 0.2e1 * t68 * pkin(3); -t75 + (MDP(17) + t77) * t43; t42 * t77 - t72 + t46 * MDP(17) - t47 * MDP(18) + (t64 * t69 - t73) * t62; t53 * t77 + (-MDP(18) * pkin(8) + MDP(15)) * t64 + (-MDP(17) * pkin(8) + t69) * t61; MDP(20) * pkin(4) ^ 2 + MDP(16); t48 * MDP(20); t51 * MDP(20); t74; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
