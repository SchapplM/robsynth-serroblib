% Calculate joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR3_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:37
% EndTime: 2022-01-23 09:14:37
% DurationCPUTime: 0.13s
% Computational Cost: add. (187->47), mult. (329->74), div. (0->0), fcn. (351->8), ass. (0->33)
t62 = sin(pkin(9));
t64 = cos(pkin(9));
t67 = sin(qJ(4));
t79 = cos(qJ(4));
t54 = t79 * t62 + t67 * t64;
t53 = t62 * t67 - t79 * t64;
t75 = t53 * MDP(13);
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t46 = t68 * t53 + t54 * t66;
t42 = t46 * MDP(20);
t47 = -t53 * t66 + t54 * t68;
t78 = -t47 * MDP(21) - t42;
t83 = -t54 * MDP(14) - t75 + t78;
t65 = cos(pkin(8));
t59 = -pkin(1) * t65 - pkin(2);
t55 = -pkin(3) * t64 + t59;
t82 = 0.2e1 * pkin(4) * t53 + 0.2e1 * t55;
t81 = 0.2e1 * t55;
t63 = sin(pkin(8));
t57 = pkin(1) * t63 + qJ(3);
t80 = pkin(6) + t57;
t77 = t62 ^ 2 + t64 ^ 2;
t76 = t64 * MDP(5);
t51 = t80 * t62;
t52 = t80 * t64;
t74 = -t79 * t51 - t52 * t67;
t40 = -pkin(7) * t54 + t74;
t71 = t67 * t51 - t79 * t52;
t41 = -t53 * pkin(7) - t71;
t73 = t47 * MDP(17) - t46 * MDP(18) + (t40 * t68 - t41 * t66) * MDP(20) + (-t40 * t66 - t41 * t68) * MDP(21);
t70 = (MDP(20) * t68 - MDP(21) * t66) * pkin(4);
t1 = [MDP(1) - 0.2e1 * t59 * t76 + (t77 * t57 ^ 2 + t59 ^ 2) * MDP(7) + t75 * t81 + t42 * t82 + (t63 ^ 2 + t65 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t77 * MDP(6) * t57 + (MDP(14) * t81 + MDP(8) * t54 - 0.2e1 * t53 * MDP(9)) * t54 + (MDP(15) * t47 - 0.2e1 * t46 * MDP(16) + MDP(21) * t82) * t47; 0; t77 * MDP(7) + MDP(4); MDP(7) * t59 - t76 - t83; 0; MDP(7); t54 * MDP(10) - t53 * MDP(11) + t74 * MDP(13) + t71 * MDP(14) + t73; t83; 0; MDP(12) + MDP(19) + 0.2e1 * t70; t73; t78; 0; MDP(19) + t70; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
