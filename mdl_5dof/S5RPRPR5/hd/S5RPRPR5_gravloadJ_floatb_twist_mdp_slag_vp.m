% Calculate Gravitation load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:02
% EndTime: 2020-01-03 11:43:04
% DurationCPUTime: 0.21s
% Computational Cost: add. (139->52), mult. (182->77), div. (0->0), fcn. (173->8), ass. (0->33)
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t70 = cos(qJ(3));
t89 = -t65 * (-qJ(4) - pkin(6)) + (t70 * pkin(3) + pkin(2)) * t66;
t87 = g(1) * t65;
t61 = qJ(3) + pkin(9) + qJ(5);
t58 = sin(t61);
t69 = sin(qJ(1));
t82 = t69 * t58;
t59 = cos(t61);
t81 = t69 * t59;
t68 = sin(qJ(3));
t80 = t69 * t68;
t79 = t69 * t70;
t71 = cos(qJ(1));
t78 = t71 * t58;
t77 = t71 * t59;
t76 = t71 * t68;
t75 = t71 * t70;
t47 = -t66 * t82 - t77;
t48 = t66 * t81 - t78;
t49 = t66 * t78 - t81;
t50 = t66 * t77 + t82;
t74 = (-g(2) * t47 - g(3) * t49 + t58 * t87) * MDP(22) + (g(2) * t48 - g(3) * t50 + t59 * t87) * MDP(23);
t73 = t71 * pkin(1) + t69 * qJ(2);
t57 = g(2) * t71 + g(3) * t69;
t56 = g(2) * t69 - g(3) * t71;
t53 = t66 * t76 - t79;
t51 = -t66 * t80 - t75;
t63 = t69 * pkin(1);
t54 = t66 * t75 + t80;
t52 = t66 * t79 - t76;
t1 = [(-g(2) * t73 - g(3) * (-t71 * qJ(2) + t63)) * MDP(7) + (-g(2) * t54 - g(3) * t52) * MDP(13) + (g(2) * t53 - g(3) * t51) * MDP(14) + (-g(2) * (pkin(3) * t80 + t73) - g(3) * (t89 * t69 + t63) + (-g(2) * t89 - g(3) * (-pkin(3) * t68 - qJ(2))) * t71) * MDP(16) + (-g(2) * t50 - g(3) * t48) * MDP(22) + (g(2) * t49 - g(3) * t47) * MDP(23) + (MDP(3) - MDP(6)) * t56 + (-t66 * MDP(4) - MDP(2) + (MDP(5) - MDP(15)) * t65) * t57; (MDP(16) + MDP(7)) * t57; (g(2) * t52 - g(3) * t54 + t70 * t87) * MDP(14) + t74 + (pkin(3) * MDP(16) + MDP(13)) * (-g(2) * t51 - g(3) * t53 + t68 * t87); (g(1) * t66 - t56 * t65) * MDP(16); t74;];
taug = t1;
