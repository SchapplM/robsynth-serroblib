% Calculate Gravitation load on the joints for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (92->36), mult. (136->52), div. (0->0), fcn. (130->8), ass. (0->23)
t61 = qJ(3) + pkin(8);
t45 = sin(t61);
t52 = sin(qJ(1));
t55 = cos(qJ(1));
t60 = cos(t61);
t34 = t55 * t45 - t52 * t60;
t35 = t52 * t45 + t55 * t60;
t54 = cos(qJ(3));
t51 = sin(qJ(3));
t66 = t52 * t51;
t36 = -t55 * t54 - t66;
t65 = t55 * t51;
t37 = -t52 * t54 + t65;
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t69 = (t53 * MDP(16) - t50 * MDP(17)) * (g(1) * t34 + g(2) * t35) - (g(1) * t36 + g(2) * t37) * MDP(9);
t64 = t55 * pkin(1) + t52 * qJ(2);
t59 = g(1) * t35 - g(2) * t34;
t57 = g(1) * t37 - g(2) * t36;
t47 = t55 * qJ(2);
t44 = t54 * pkin(3) + pkin(2);
t38 = g(1) * t52 - g(2) * t55;
t1 = [(-g(1) * (-t52 * pkin(1) + t47) - g(2) * t64) * MDP(6) - t57 * MDP(8) + (-g(1) * (pkin(3) * t65 + t47 + (-pkin(1) - t44) * t52) - g(2) * (pkin(3) * t66 + t55 * t44 + t64)) * MDP(10) + (MDP(3) - MDP(5)) * (g(1) * t55 + g(2) * t52) + (MDP(2) + MDP(4)) * t38 - t69; (-MDP(10) - MDP(6)) * t38; (MDP(10) * pkin(3) + MDP(8)) * t57 + t69; g(3) * MDP(10); (g(3) * t53 + t59 * t50) * MDP(16) + (-g(3) * t50 + t59 * t53) * MDP(17);];
taug = t1;
