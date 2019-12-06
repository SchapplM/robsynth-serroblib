% Calculate Gravitation load on the joints for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:09
% EndTime: 2019-12-05 16:10:11
% DurationCPUTime: 0.29s
% Computational Cost: add. (140->58), mult. (232->85), div. (0->0), fcn. (212->8), ass. (0->34)
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t67 = g(1) * t56 + g(2) * t55;
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t45 = -g(3) * t61 + t67 * t59;
t79 = g(3) * t59;
t60 = cos(qJ(3));
t77 = t55 * t60;
t76 = t55 * t61;
t75 = t56 * t60;
t74 = t56 * t61;
t57 = -qJ(4) - pkin(6);
t73 = t57 * t61;
t58 = sin(qJ(3));
t72 = t58 * t61;
t71 = t60 * t61;
t70 = -MDP(13) - MDP(17);
t69 = t56 * t72;
t50 = pkin(3) * t60 + pkin(2);
t68 = t61 * t50 - t59 * t57;
t54 = qJ(3) + pkin(8);
t51 = sin(t54);
t52 = cos(t54);
t66 = pkin(4) * t52 + qJ(5) * t51;
t64 = -t55 * t72 - t75;
t41 = t51 * t76 + t56 * t52;
t43 = t51 * t74 - t55 * t52;
t63 = g(1) * t43 + g(2) * t41 + t51 * t79;
t46 = t67 * t61 + t79;
t49 = pkin(3) * t77;
t44 = t55 * t51 + t52 * t74;
t42 = -t56 * t51 + t52 * t76;
t1 = [(-MDP(1) + t70) * g(3); (-g(3) * t68 + t67 * (t50 * t59 + t73)) * MDP(13) + (-g(3) * (t66 * t61 + t68) + t67 * (t73 - (-t50 - t66) * t59)) * MDP(17) + (MDP(4) - MDP(12) - MDP(15)) * t46 + (t60 * MDP(10) - t58 * MDP(11) + t52 * MDP(14) + t51 * MDP(16) + MDP(3)) * t45; (-g(1) * (-t69 + t77) - g(2) * t64 + t58 * t79) * MDP(10) + (-g(1) * (-t55 * t58 - t56 * t71) - g(2) * (-t55 * t71 + t56 * t58) + t60 * t79) * MDP(11) + (-g(1) * t49 + (g(2) * t75 + t46 * t58) * pkin(3)) * MDP(13) + t63 * MDP(14) + (-g(1) * t44 - g(2) * t42 - t52 * t79) * MDP(16) + (-g(1) * (-pkin(3) * t69 - t43 * pkin(4) + t44 * qJ(5) + t49) - g(2) * (t64 * pkin(3) - t41 * pkin(4) + t42 * qJ(5)) - (-pkin(3) * t58 - pkin(4) * t51 + qJ(5) * t52) * t79) * MDP(17); t70 * t45; -t63 * MDP(17);];
taug = t1;
