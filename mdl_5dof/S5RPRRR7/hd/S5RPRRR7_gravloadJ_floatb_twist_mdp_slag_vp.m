% Calculate Gravitation load on the joints for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:09
% EndTime: 2019-12-31 19:04:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (161->45), mult. (172->73), div. (0->0), fcn. (168->10), ass. (0->29)
t51 = sin(qJ(3));
t70 = t51 * MDP(11);
t49 = qJ(4) + qJ(5);
t46 = sin(t49);
t47 = cos(t49);
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t69 = t53 * MDP(17) - t50 * MDP(18) + t47 * MDP(24) - t46 * MDP(25) + MDP(10);
t68 = g(3) * t51;
t54 = cos(qJ(3));
t67 = t46 * t54;
t66 = t47 * t54;
t65 = t50 * t54;
t64 = t53 * t54;
t48 = qJ(1) + pkin(9);
t44 = sin(t48);
t45 = cos(t48);
t36 = t44 * t67 + t45 * t47;
t37 = -t44 * t66 + t45 * t46;
t38 = t44 * t47 - t45 * t67;
t39 = t44 * t46 + t45 * t66;
t63 = (-g(1) * t38 + g(2) * t36 + t46 * t68) * MDP(24) + (g(1) * t39 - g(2) * t37 + t47 * t68) * MDP(25);
t55 = cos(qJ(1));
t52 = sin(qJ(1));
t43 = t44 * t50 + t45 * t64;
t42 = t44 * t53 - t45 * t65;
t41 = -t44 * t64 + t45 * t50;
t40 = t44 * t65 + t45 * t53;
t1 = [(g(1) * t55 + g(2) * t52) * MDP(3) + (-g(1) * t41 - g(2) * t43) * MDP(17) + (-g(1) * t40 - g(2) * t42) * MDP(18) + (-g(1) * t37 - g(2) * t39) * MDP(24) + (-g(1) * t36 - g(2) * t38) * MDP(25) + (t54 * MDP(10) - t70) * (g(1) * t44 - g(2) * t45) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t52 - g(2) * t55); -g(3) * MDP(4); (-t69 * t54 + t70) * g(3) + (MDP(11) * t54 + t69 * t51) * (g(1) * t45 + g(2) * t44); (-g(1) * t42 + g(2) * t40 + t50 * t68) * MDP(17) + (g(1) * t43 - g(2) * t41 + t53 * t68) * MDP(18) + t63; t63;];
taug = t1;
