% Calculate joint inertia matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRR5_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:39
% EndTime: 2019-12-31 17:56:39
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->37), mult. (133->46), div. (0->0), fcn. (90->6), ass. (0->25)
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t37 = sin(qJ(5));
t39 = cos(qJ(5));
t48 = t39 * MDP(16);
t45 = -t37 * MDP(17) + t48;
t55 = (MDP(9) + t45) * t40 - t38 * MDP(10);
t36 = cos(pkin(8));
t31 = t36 * pkin(1) + pkin(2);
t29 = -pkin(3) - t31;
t35 = sin(pkin(8));
t30 = t35 * pkin(1) + qJ(3);
t26 = -t40 * t29 + t38 * t30;
t24 = pkin(4) + t26;
t54 = pkin(4) + t24;
t53 = t26 * MDP(9);
t52 = t37 ^ 2 * MDP(11) + MDP(8);
t27 = t38 * t29 + t40 * t30;
t51 = t27 * MDP(10);
t49 = t39 * MDP(12);
t47 = 0.2e1 * t37 * t49 + t52;
t46 = -t37 * MDP(13) - t39 * MDP(14);
t44 = -MDP(16) * t37 - MDP(17) * t39;
t42 = 0.2e1 * t45;
t1 = [MDP(1) + (t30 ^ 2 + t31 ^ 2) * MDP(7) + (t35 ^ 2 + t36 ^ 2) * MDP(4) * pkin(1) ^ 2 + t24 * t42 + 0.2e1 * t51 + 0.2e1 * t31 * MDP(5) + 0.2e1 * t30 * MDP(6) + 0.2e1 * t53 + t47; 0; MDP(4) + MDP(7); -t31 * MDP(7) - MDP(5) - t55; 0; MDP(7); -t51 - t53 - t54 * t48 + (t54 * MDP(17) - 0.2e1 * t49) * t37 - t52; 0; t55; pkin(4) * t42 + t47; t44 * (-pkin(7) + t27) + t46; -t45; t44 * t38; t44 * pkin(7) - t46; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
