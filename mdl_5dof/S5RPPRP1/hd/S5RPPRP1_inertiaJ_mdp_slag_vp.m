% Calculate joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP1_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:46
% EndTime: 2022-01-23 09:12:47
% DurationCPUTime: 0.21s
% Computational Cost: add. (189->69), mult. (342->98), div. (0->0), fcn. (272->6), ass. (0->31)
t54 = sin(qJ(4));
t62 = MDP(14) + MDP(16);
t69 = t62 * t54;
t51 = sin(pkin(7));
t44 = t51 * pkin(1) + qJ(3);
t68 = t44 * t54;
t55 = cos(qJ(4));
t49 = t55 ^ 2;
t67 = t54 ^ 2 + t49;
t66 = MDP(18) * pkin(4);
t50 = sin(pkin(8));
t65 = qJ(5) * t50;
t64 = t54 * MDP(11);
t63 = -MDP(13) - MDP(15);
t52 = cos(pkin(8));
t61 = t55 * t52 * t44;
t53 = cos(pkin(7));
t45 = -t53 * pkin(1) - pkin(2);
t42 = -t52 * pkin(3) - t50 * pkin(6) + t45;
t40 = t55 * t42;
t36 = -t55 * t65 + t40 + (-pkin(4) - t68) * t52;
t37 = t61 + (t42 - t65) * t54;
t60 = t36 * t55 + t37 * t54;
t59 = -t63 + t66;
t58 = t54 * MDP(15) + t55 * MDP(16);
t57 = (-t52 * t68 + t40) * MDP(13) - (t54 * t42 + t61) * MDP(14) - t37 * MDP(16);
t47 = t52 ^ 2;
t46 = t50 ^ 2;
t43 = t44 ^ 2;
t41 = (pkin(4) * t54 + t44) * t50;
t1 = [MDP(1) + (t47 * t43 + t45 ^ 2) * MDP(7) + t47 * MDP(12) + (t36 ^ 2 + t37 ^ 2 + t41 ^ 2) * MDP(18) + (t51 ^ 2 + t53 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t55 * t54 * MDP(9) + t43 * MDP(7) + t49 * MDP(8)) * t46 + 0.2e1 * (-t60 * MDP(17) + t58 * t41) * t50 + 0.2e1 * (t47 * MDP(6) + (t54 * MDP(13) + t55 * MDP(14) + MDP(6)) * t46) * t44 + 0.2e1 * (-t45 * MDP(5) + (-t55 * MDP(10) + t64) * t50 - t36 * MDP(15) - t57) * t52; (-t41 * t52 + (-t36 * t54 + t37 * t55) * t50) * MDP(18); MDP(4) + (t46 + t47) * MDP(7) + (t67 * t46 + t47) * MDP(18); t45 * MDP(7) + t60 * MDP(18) - t67 * MDP(17) * t50 + (t63 * t55 - MDP(5) + t69) * t52; 0; t67 * MDP(18) + MDP(7); t36 * t66 + t40 * MDP(15) + (-MDP(12) + (-0.2e1 * pkin(4) - t68) * MDP(15)) * t52 + (-t64 + (-MDP(15) * qJ(5) - MDP(17) * pkin(4) + MDP(10)) * t55) * t50 + t57; (-t59 * t54 - t62 * t55) * t50; t59 * t55 - t69; MDP(12) + (0.2e1 * MDP(15) + t66) * pkin(4); t41 * MDP(18) + t58 * t50; -t52 * MDP(18); 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
