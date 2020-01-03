% Calculate Gravitation load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:37
% EndTime: 2020-01-03 12:07:37
% DurationCPUTime: 0.06s
% Computational Cost: add. (147->31), mult. (100->46), div. (0->0), fcn. (70->10), ass. (0->20)
t46 = qJ(1) + qJ(2);
t45 = qJ(3) + t46;
t41 = sin(t45);
t43 = sin(t46);
t57 = pkin(2) * t43 + pkin(3) * t41;
t42 = cos(t45);
t44 = cos(t46);
t56 = pkin(2) * t44 + pkin(3) * t42;
t40 = pkin(9) + t45;
t34 = sin(t40);
t35 = cos(t40);
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t52 = -g(2) * t42 - g(3) * t41;
t55 = (g(2) * t41 - g(3) * t42) * MDP(9) + t52 * MDP(8) + (-t49 * MDP(16) + t47 * MDP(17)) * (g(2) * t35 + g(3) * t34);
t53 = g(2) * t34 - g(3) * t35;
t51 = (g(2) * t43 - g(3) * t44) * MDP(6) + (-g(2) * t44 - g(3) * t43) * MDP(5) + t55;
t50 = cos(qJ(1));
t48 = sin(qJ(1));
t1 = [(-g(2) * t50 - g(3) * t48) * MDP(2) + (g(2) * t48 - g(3) * t50) * MDP(3) + (-g(2) * (t50 * pkin(1) + t56) - g(3) * (t48 * pkin(1) + t57)) * MDP(10) + t51; (-g(2) * t56 - g(3) * t57) * MDP(10) + t51; t52 * MDP(10) * pkin(3) + t55; -g(1) * MDP(10); (-g(1) * t49 + t53 * t47) * MDP(16) + (g(1) * t47 + t53 * t49) * MDP(17);];
taug = t1;
