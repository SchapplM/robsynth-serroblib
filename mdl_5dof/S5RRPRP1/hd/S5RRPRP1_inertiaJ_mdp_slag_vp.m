% Calculate joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:20:00
% EndTime: 2022-01-20 10:20:01
% DurationCPUTime: 0.15s
% Computational Cost: add. (230->76), mult. (364->99), div. (0->0), fcn. (278->6), ass. (0->43)
t64 = sin(qJ(4));
t57 = t64 * MDP(16);
t66 = cos(qJ(4));
t73 = -t66 * MDP(15) + t57;
t88 = 2 * MDP(17);
t65 = sin(qJ(2));
t87 = pkin(1) * t65;
t86 = t66 * pkin(4);
t67 = cos(qJ(2));
t56 = t67 * pkin(1) + pkin(2);
t62 = sin(pkin(8));
t63 = cos(pkin(8));
t45 = t63 * t56 - t62 * t87;
t42 = -pkin(3) - t45;
t41 = t42 - t86;
t55 = -t63 * pkin(2) - pkin(3);
t49 = t55 - t86;
t85 = t41 + t49;
t84 = t42 + t55;
t46 = t62 * t56 + t63 * t87;
t83 = t64 * MDP(10) + t66 * MDP(11);
t82 = MDP(18) * pkin(4);
t43 = pkin(7) + t46;
t81 = qJ(5) + t43;
t54 = t62 * pkin(2) + pkin(7);
t80 = qJ(5) + t54;
t79 = t41 * MDP(18);
t78 = t49 * MDP(18);
t77 = t64 * MDP(14);
t76 = t64 * MDP(17);
t61 = t64 ^ 2;
t74 = 0.2e1 * t64 * t66 * MDP(9) + t61 * MDP(8) + MDP(4);
t72 = -MDP(13) * t64 - MDP(14) * t66;
t71 = 0.2e1 * t73;
t70 = (t67 * MDP(5) - t65 * MDP(6)) * pkin(1);
t69 = -0.2e1 * t66 * MDP(13) + 0.2e1 * t77;
t48 = t80 * t66;
t47 = t80 * t64;
t44 = t48 * t66;
t40 = t81 * t66;
t39 = t81 * t64;
t38 = t40 * t66;
t1 = [MDP(1) + (t45 ^ 2 + t46 ^ 2) * MDP(7) + (t39 * t64 + t38) * t88 + (t39 ^ 2 + t40 ^ 2) * MDP(18) + t42 * t69 + (t71 + t79) * t41 + 0.2e1 * t70 + t74; (t38 + t44) * MDP(17) + (t39 * t47 + t40 * t48 + t41 * t49) * MDP(18) + (t45 * t63 + t46 * t62) * MDP(7) * pkin(2) + t70 + (-t84 * MDP(13) - t85 * MDP(15)) * t66 + (t84 * MDP(14) + t85 * MDP(16) + (t39 + t47) * MDP(17)) * t64 + t74; (t47 * t64 + t44) * t88 + (t47 ^ 2 + t48 ^ 2) * MDP(18) + (t62 ^ 2 + t63 ^ 2) * MDP(7) * pkin(2) ^ 2 + t55 * t69 + (t71 + t78) * t49 + t74; (-t39 * t66 + t40 * t64) * MDP(18); (-t47 * t66 + t48 * t64) * MDP(18); MDP(7) + (t66 ^ 2 + t61) * MDP(18); -t39 * MDP(15) - t40 * MDP(16) + t72 * t43 + (-t39 * MDP(18) - t76) * pkin(4) + t83; -t47 * MDP(15) - t48 * MDP(16) + t72 * t54 + (-MDP(18) * t47 - t76) * pkin(4) + t83; -t77 - t57 + (MDP(13) + MDP(15) + t82) * t66; MDP(12) + (0.2e1 * MDP(15) + t82) * pkin(4); t73 + t79; t73 + t78; 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
