% Calculate Gravitation load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:02
% EndTime: 2021-01-15 11:24:03
% DurationCPUTime: 0.21s
% Computational Cost: add. (123->51), mult. (174->64), div. (0->0), fcn. (135->8), ass. (0->25)
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t67 = -g(1) * t53 + g(2) * t55;
t49 = qJ(3) + pkin(7);
t44 = sin(t49);
t64 = g(3) * t44;
t52 = sin(qJ(3));
t63 = g(3) * t52;
t48 = pkin(1) + pkin(6) + qJ(4);
t62 = t48 * t53;
t61 = MDP(14) + MDP(18);
t60 = MDP(15) - MDP(20);
t59 = -MDP(17) - MDP(21);
t41 = g(1) * t55 + g(2) * t53;
t50 = sin(pkin(7));
t51 = cos(pkin(7));
t38 = pkin(4) * t51 + qJ(5) * t50 + pkin(3);
t39 = -t50 * pkin(4) + qJ(5) * t51;
t54 = cos(qJ(3));
t57 = t38 * t52 - t39 * t54 + qJ(2);
t56 = t54 * t67 + t63;
t45 = cos(t49);
t43 = pkin(3) * t52 + qJ(2);
t42 = t48 * t55;
t1 = [(-g(1) * (-t53 * pkin(1) + qJ(2) * t55) - g(2) * (pkin(1) * t55 + t53 * qJ(2))) * MDP(6) + (-g(1) * (t43 * t55 - t62) - g(2) * (t43 * t53 + t42)) * MDP(17) + (-g(1) * (t57 * t55 - t62) - g(2) * (t57 * t53 + t42)) * MDP(21) - (MDP(2) - MDP(4) + MDP(16) + MDP(19)) * t67 + (-t52 * MDP(12) - t54 * MDP(13) - t61 * t44 - t60 * t45 + MDP(3) - MDP(5)) * t41; -(-MDP(6) + t59) * t67; t56 * MDP(12) + (g(3) * t54 - t52 * t67) * MDP(13) + (-g(3) * (-pkin(4) * t44 + qJ(5) * t45) + t67 * (t38 * t54 + t39 * t52)) * MDP(21) + t61 * (t45 * t67 + t64) + t60 * (g(3) * t45 - t44 * t67) + (t56 * MDP(17) + MDP(21) * t63) * pkin(3); t59 * t41; (-t64 - t67 * (-t50 * t52 + t51 * t54)) * MDP(21);];
taug = t1;
