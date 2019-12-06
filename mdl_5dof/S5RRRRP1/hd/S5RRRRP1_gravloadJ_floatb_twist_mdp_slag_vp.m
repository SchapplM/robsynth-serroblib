% Calculate Gravitation load on the joints for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (175->39), mult. (163->52), div. (0->0), fcn. (126->8), ass. (0->21)
t45 = qJ(2) + qJ(3);
t43 = qJ(4) + t45;
t38 = sin(t43);
t39 = cos(t43);
t47 = sin(qJ(1));
t49 = cos(qJ(1));
t51 = g(1) * t49 + g(2) * t47;
t50 = -g(3) * t39 + t51 * t38;
t56 = t50 * MDP(23) + (g(3) * t38 + t51 * t39) * MDP(24);
t41 = cos(t45);
t55 = pkin(3) * t41 + pkin(4) * t39;
t48 = cos(qJ(2));
t54 = t48 * pkin(2) + t55;
t40 = sin(t45);
t53 = (-g(3) * t41 + t51 * t40) * MDP(16) + (g(3) * t40 + t51 * t41) * MDP(17) + t56;
t52 = -pkin(3) * t40 - pkin(4) * t38;
t34 = g(1) * t47 - g(2) * t49;
t46 = sin(qJ(2));
t42 = -qJ(5) - pkin(8) - pkin(7) - pkin(6);
t31 = pkin(1) + t54;
t1 = [(-g(1) * (-t47 * t31 - t49 * t42) - g(2) * (t49 * t31 - t47 * t42)) * MDP(26) + (MDP(3) - MDP(25)) * t51 + (-t46 * MDP(10) + MDP(16) * t41 - MDP(17) * t40 + MDP(23) * t39 - MDP(24) * t38 + t48 * MDP(9) + MDP(2)) * t34; (-g(3) * t48 + t51 * t46) * MDP(9) + (g(3) * t46 + t51 * t48) * MDP(10) + (-g(3) * t54 - t51 * (-t46 * pkin(2) + t52)) * MDP(26) + t53; (-g(3) * t55 - t51 * t52) * MDP(26) + t53; t50 * MDP(26) * pkin(4) + t56; -t34 * MDP(26);];
taug = t1;
