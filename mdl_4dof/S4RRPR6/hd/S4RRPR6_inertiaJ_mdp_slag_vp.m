% Calculate joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR6_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:23
% EndTime: 2021-01-15 10:46:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (171->57), mult. (327->90), div. (0->0), fcn. (325->6), ass. (0->30)
t56 = sin(pkin(7));
t57 = cos(pkin(7));
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t48 = t56 * t59 - t57 * t61;
t55 = -t61 * pkin(2) - pkin(1);
t73 = 0.2e1 * t48 * pkin(3) + 0.2e1 * t55;
t72 = 0.2e1 * t61;
t71 = pkin(2) * t56;
t70 = -qJ(3) - pkin(5);
t49 = t56 * t61 + t57 * t59;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t40 = t60 * t48 + t58 * t49;
t69 = t40 * MDP(20);
t54 = t57 * pkin(2) + pkin(3);
t68 = (t60 * t54 - t58 * t71) * MDP(20);
t67 = (-t58 * t54 - t60 * t71) * MDP(21);
t66 = t48 * MDP(11);
t65 = t49 * MDP(12);
t64 = t55 * MDP(14);
t51 = t70 * t59;
t52 = t70 * t61;
t42 = t57 * t51 + t56 * t52;
t36 = -t49 * pkin(6) + t42;
t43 = t56 * t51 - t57 * t52;
t37 = -t48 * pkin(6) + t43;
t41 = -t58 * t48 + t60 * t49;
t63 = t41 * MDP(17) - t40 * MDP(18) + (t60 * t36 - t58 * t37) * MDP(20) + (-t58 * t36 - t60 * t37) * MDP(21);
t1 = [MDP(1) + pkin(1) * MDP(9) * t72 + 0.2e1 * (-t42 * t49 - t43 * t48) * MDP(13) + (t42 ^ 2 + t43 ^ 2) * MDP(14) + t69 * t73 + (t64 + 0.2e1 * t65 + 0.2e1 * t66) * t55 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t59 + MDP(5) * t72) * t59 + (MDP(15) * t41 - 0.2e1 * t40 * MDP(16) + MDP(21) * t73) * t41; t42 * MDP(11) - t43 * MDP(12) + t59 * MDP(6) + t61 * MDP(7) + (-t61 * MDP(10) - t59 * MDP(9)) * pkin(5) + ((-t48 * t56 - t49 * t57) * MDP(13) + (t42 * t57 + t43 * t56) * MDP(14)) * pkin(2) + t63; MDP(8) + MDP(19) + 0.2e1 * t68 + 0.2e1 * t67 + (0.2e1 * t57 * MDP(11) - 0.2e1 * t56 * MDP(12) + (t56 ^ 2 + t57 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t41 * MDP(21) + t64 + t65 + t66 + t69; 0; MDP(14); t63; MDP(19) + t67 + t68; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
