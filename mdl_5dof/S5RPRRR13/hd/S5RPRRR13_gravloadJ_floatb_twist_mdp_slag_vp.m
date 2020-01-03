% Calculate Gravitation load on the joints for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:28
% EndTime: 2019-12-31 19:15:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (111->47), mult. (182->70), div. (0->0), fcn. (176->8), ass. (0->31)
t55 = cos(qJ(3));
t73 = MDP(13) * t55;
t50 = qJ(4) + qJ(5);
t47 = sin(t50);
t48 = cos(t50);
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t72 = t54 * MDP(19) - t51 * MDP(20) + t48 * MDP(26) - t47 * MDP(27) + MDP(12);
t71 = g(3) * t55;
t53 = sin(qJ(1));
t70 = t53 * t47;
t69 = t53 * t48;
t68 = t53 * t51;
t67 = t53 * t54;
t56 = cos(qJ(1));
t66 = t56 * t47;
t65 = t56 * t48;
t64 = t56 * t51;
t63 = t56 * t54;
t52 = sin(qJ(3));
t37 = -t52 * t70 + t65;
t38 = t52 * t69 + t66;
t39 = t52 * t66 + t69;
t40 = t52 * t65 - t70;
t62 = (-g(1) * t37 - g(2) * t39 + t47 * t71) * MDP(26) + (g(1) * t38 - g(2) * t40 + t48 * t71) * MDP(27);
t45 = g(1) * t53 - g(2) * t56;
t44 = t52 * t63 - t68;
t43 = t52 * t64 + t67;
t42 = t52 * t67 + t64;
t41 = -t52 * t68 + t63;
t1 = [(-g(1) * (-t53 * pkin(1) + t56 * qJ(2)) - g(2) * (t56 * pkin(1) + t53 * qJ(2))) * MDP(6) + (-g(1) * t44 - g(2) * t42) * MDP(19) + (g(1) * t43 - g(2) * t41) * MDP(20) + (-g(1) * t40 - g(2) * t38) * MDP(26) + (g(1) * t39 - g(2) * t37) * MDP(27) + (MDP(2) - MDP(4)) * t45 + (-MDP(12) * t52 + MDP(3) - MDP(5) - t73) * (g(1) * t56 + g(2) * t53); -t45 * MDP(6); (t72 * t52 + t73) * g(3) + (MDP(13) * t52 - t72 * t55) * t45; (-g(1) * t41 - g(2) * t43 + t51 * t71) * MDP(19) + (g(1) * t42 - g(2) * t44 + t54 * t71) * MDP(20) + t62; t62;];
taug = t1;
