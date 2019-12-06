% Calculate Gravitation load on the joints for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:38
% EndTime: 2019-12-05 14:59:38
% DurationCPUTime: 0.08s
% Computational Cost: add. (63->33), mult. (161->63), div. (0->0), fcn. (189->10), ass. (0->25)
t40 = sin(pkin(9));
t41 = sin(pkin(8));
t57 = t40 * t41;
t47 = sin(qJ(4));
t56 = t41 * t47;
t49 = cos(qJ(4));
t55 = t41 * t49;
t42 = sin(pkin(7));
t44 = cos(pkin(8));
t54 = t42 * t44;
t45 = cos(pkin(7));
t53 = t45 * t40;
t43 = cos(pkin(9));
t52 = t45 * t43;
t51 = MDP(2) + MDP(3);
t48 = cos(qJ(5));
t46 = sin(qJ(5));
t38 = t43 * t55 - t44 * t47;
t36 = t42 * t40 + t44 * t52;
t35 = -t42 * t43 + t44 * t53;
t34 = t43 * t54 - t53;
t33 = t40 * t54 + t52;
t32 = t36 * t49 + t45 * t56;
t30 = t34 * t49 + t42 * t56;
t1 = [(-MDP(1) - t51) * g(3); t51 * (-g(1) * t42 + g(2) * t45); (g(3) * t44 + (-g(1) * t45 - g(2) * t42) * t41) * MDP(3); (g(1) * t32 + g(2) * t30 + g(3) * t38) * MDP(6) + (-MDP(12) * t48 + MDP(13) * t46 - MDP(5)) * (g(1) * (-t36 * t47 + t45 * t55) + g(2) * (-t34 * t47 + t42 * t55) + g(3) * (-t43 * t56 - t44 * t49)); (-g(1) * (-t32 * t46 + t35 * t48) - g(2) * (-t30 * t46 + t33 * t48) - g(3) * (-t38 * t46 + t48 * t57)) * MDP(12) + (-g(1) * (-t32 * t48 - t35 * t46) - g(2) * (-t30 * t48 - t33 * t46) - g(3) * (-t38 * t48 - t46 * t57)) * MDP(13);];
taug = t1;
