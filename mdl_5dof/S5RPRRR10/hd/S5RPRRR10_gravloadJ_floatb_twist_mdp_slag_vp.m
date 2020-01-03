% Calculate Gravitation load on the joints for
% S5RPRRR10
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:45
% EndTime: 2019-12-31 19:10:45
% DurationCPUTime: 0.17s
% Computational Cost: add. (161->49), mult. (189->72), div. (0->0), fcn. (182->10), ass. (0->33)
t53 = pkin(9) + qJ(3);
t49 = sin(t53);
t77 = MDP(14) * t49;
t54 = qJ(4) + qJ(5);
t51 = sin(t54);
t52 = cos(t54);
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t76 = t59 * MDP(20) - t57 * MDP(21) + t52 * MDP(27) - t51 * MDP(28) + MDP(13);
t75 = g(3) * t49;
t58 = sin(qJ(1));
t74 = t58 * t51;
t73 = t58 * t52;
t72 = t58 * t57;
t71 = t58 * t59;
t60 = cos(qJ(1));
t70 = t60 * t51;
t69 = t60 * t52;
t68 = t60 * t57;
t67 = t60 * t59;
t50 = cos(t53);
t39 = t50 * t74 + t69;
t40 = -t50 * t73 + t70;
t41 = -t50 * t70 + t73;
t42 = t50 * t69 + t74;
t66 = (-g(1) * t41 + g(2) * t39 + t51 * t75) * MDP(27) + (g(1) * t42 - g(2) * t40 + t52 * t75) * MDP(28);
t61 = g(1) * t60 + g(2) * t58;
t47 = g(1) * t58 - g(2) * t60;
t46 = t50 * t67 + t72;
t45 = -t50 * t68 + t71;
t44 = -t50 * t71 + t68;
t43 = t50 * t72 + t67;
t1 = [(-g(1) * (-t58 * pkin(1) + t60 * qJ(2)) - g(2) * (t60 * pkin(1) + t58 * qJ(2))) * MDP(7) + (-g(1) * t44 - g(2) * t46) * MDP(20) + (-g(1) * t43 - g(2) * t45) * MDP(21) + (-g(1) * t40 - g(2) * t42) * MDP(27) + (-g(1) * t39 - g(2) * t41) * MDP(28) + (MDP(3) - MDP(6)) * t61 + (MDP(13) * t50 - t77 + MDP(4) * cos(pkin(9)) - MDP(5) * sin(pkin(9)) + MDP(2)) * t47; -t47 * MDP(7); (-t76 * t50 + t77) * g(3) + (MDP(14) * t50 + t76 * t49) * t61; (-g(1) * t45 + g(2) * t43 + t57 * t75) * MDP(20) + (g(1) * t46 - g(2) * t44 + t59 * t75) * MDP(21) + t66; t66;];
taug = t1;
