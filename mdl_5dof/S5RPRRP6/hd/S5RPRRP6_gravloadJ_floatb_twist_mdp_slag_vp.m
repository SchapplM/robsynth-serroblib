% Calculate Gravitation load on the joints for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (130->44), mult. (166->66), div. (0->0), fcn. (143->8), ass. (0->27)
t47 = qJ(1) + pkin(8);
t45 = sin(t47);
t46 = cos(t47);
t61 = g(1) * t46 + g(2) * t45;
t73 = MDP(11) - MDP(19);
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t37 = -g(3) * t53 + t61 * t50;
t67 = g(3) * t50;
t49 = sin(qJ(4));
t65 = t49 * t53;
t52 = cos(qJ(4));
t64 = t52 * t53;
t62 = pkin(4) * t49 + pkin(6);
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t59 = g(1) * t51 - g(2) * t54;
t44 = t52 * pkin(4) + pkin(3);
t48 = -qJ(5) - pkin(7);
t58 = t53 * t44 - t50 * t48;
t56 = pkin(2) + t58;
t41 = t45 * t52 - t46 * t65;
t39 = t45 * t65 + t46 * t52;
t55 = t59 * pkin(1);
t42 = t45 * t49 + t46 * t64;
t40 = -t45 * t64 + t46 * t49;
t1 = [t59 * MDP(2) + (g(1) * t54 + g(2) * t51) * MDP(3) + MDP(4) * t55 + (-g(1) * t40 - g(2) * t42) * MDP(17) + (-g(1) * t39 - g(2) * t41) * MDP(18) + (t55 + (-g(1) * t62 - g(2) * t56) * t46 + (g(1) * t56 - g(2) * t62) * t45) * MDP(20) + (t53 * MDP(10) - t73 * t50) * (g(1) * t45 - g(2) * t46); (-MDP(20) - MDP(4)) * g(3); (-g(3) * t58 + t61 * (t44 * t50 + t48 * t53)) * MDP(20) + t73 * (t61 * t53 + t67) + (MDP(17) * t52 - MDP(18) * t49 + MDP(10)) * t37; (g(1) * t42 - g(2) * t40 + t52 * t67) * MDP(18) + (pkin(4) * MDP(20) + MDP(17)) * (-g(1) * t41 + g(2) * t39 + t49 * t67); -t37 * MDP(20);];
taug = t1;
