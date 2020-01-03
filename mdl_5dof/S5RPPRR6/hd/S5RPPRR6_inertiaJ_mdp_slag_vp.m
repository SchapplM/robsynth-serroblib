% Calculate joint inertia matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR6_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:08
% EndTime: 2019-12-31 17:58:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (195->57), mult. (371->87), div. (0->0), fcn. (368->8), ass. (0->41)
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t63 = sin(qJ(4));
t81 = cos(qJ(4));
t47 = t63 * t58 - t81 * t60;
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t69 = t64 * MDP(21) - t62 * MDP(22);
t67 = MDP(14) + t69;
t48 = t81 * t58 + t63 * t60;
t73 = t48 * MDP(15);
t84 = -t67 * t47 - t73;
t61 = cos(pkin(8));
t53 = -t61 * pkin(1) - pkin(2);
t49 = -t60 * pkin(3) + t53;
t83 = 0.2e1 * t49;
t59 = sin(pkin(8));
t51 = t59 * pkin(1) + qJ(3);
t82 = pkin(6) + t51;
t44 = t82 * t58;
t45 = t82 * t60;
t42 = t81 * t44 + t63 * t45;
t80 = t42 * t48;
t79 = t62 * t64;
t78 = t58 ^ 2 + t60 ^ 2;
t77 = t53 * MDP(8);
t76 = t58 * MDP(6);
t75 = t60 * MDP(5);
t74 = t47 * MDP(20);
t72 = MDP(17) * t79;
t71 = t78 * MDP(8);
t70 = MDP(18) * t64 - MDP(19) * t62;
t68 = -MDP(21) * t62 - MDP(22) * t64;
t66 = t62 * MDP(18) + t64 * MDP(19) + t68 * pkin(7);
t57 = t64 ^ 2;
t56 = t62 ^ 2;
t43 = -t63 * t44 + t81 * t45;
t41 = t47 * pkin(4) - t48 * pkin(7) + t49;
t40 = t62 * t41 + t64 * t43;
t39 = t64 * t41 - t62 * t43;
t1 = [t73 * t83 + MDP(1) + (t59 ^ 2 + t61 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t75 + 0.2e1 * t76 + t77) * t53 + (t57 * MDP(16) + MDP(9) - 0.2e1 * t72) * t48 ^ 2 + (MDP(14) * t83 + t74 + 0.2e1 * (-MDP(10) + t70) * t48) * t47 + 0.2e1 * (t39 * t47 + t62 * t80) * MDP(21) + 0.2e1 * (-t40 * t47 + t64 * t80) * MDP(22) + (0.2e1 * t78 * MDP(7) + t71 * t51) * t51; 0; MDP(4) + t71; -t75 + t76 + t77 - t84; 0; MDP(8); -t43 * MDP(15) - t67 * t42 + (-MDP(12) + t66) * t47 + (MDP(11) + MDP(16) * t79 + (-t56 + t57) * MDP(17) + t68 * pkin(4)) * t48; t84; 0; t56 * MDP(16) + 0.2e1 * pkin(4) * t69 + MDP(13) + 0.2e1 * t72; t39 * MDP(21) - t40 * MDP(22) + t70 * t48 + t74; t68 * t48; t69; t66; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
