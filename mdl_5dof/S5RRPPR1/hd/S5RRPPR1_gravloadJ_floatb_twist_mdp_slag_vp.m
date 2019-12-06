% Calculate Gravitation load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:27
% EndTime: 2019-12-05 18:18:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (155->37), mult. (113->49), div. (0->0), fcn. (80->10), ass. (0->22)
t45 = qJ(1) + qJ(2);
t42 = sin(t45);
t59 = pkin(2) * t42;
t43 = cos(t45);
t58 = pkin(2) * t43;
t48 = sin(qJ(1));
t57 = t48 * pkin(1);
t49 = cos(qJ(1));
t56 = t49 * pkin(1);
t41 = pkin(8) + t45;
t37 = sin(t41);
t38 = cos(t41);
t55 = g(2) * t38 + g(3) * t37;
t54 = g(2) * t37 - g(3) * t38;
t53 = g(2) * t43 + g(3) * t42;
t52 = -t37 * pkin(3) + t38 * qJ(4) - t59;
t44 = pkin(9) + qJ(5);
t39 = sin(t44);
t40 = cos(t44);
t51 = t54 * MDP(10) + (-g(2) * t42 + g(3) * t43) * MDP(6) + t53 * MDP(5) + (t40 * MDP(17) - t39 * MDP(18) + cos(pkin(9)) * MDP(8) - sin(pkin(9)) * MDP(9)) * t55;
t50 = -t38 * pkin(3) - t37 * qJ(4) - t58;
t1 = [(g(2) * t49 + g(3) * t48) * MDP(2) + (-g(2) * t48 + g(3) * t49) * MDP(3) + (-g(2) * (-t56 - t58) - g(3) * (-t57 - t59)) * MDP(7) + (-g(2) * (t50 - t56) - g(3) * (t52 - t57)) * MDP(11) + t51; t53 * pkin(2) * MDP(7) + (-g(2) * t50 - g(3) * t52) * MDP(11) + t51; (-MDP(11) - MDP(7)) * g(1); -t55 * MDP(11); (-g(1) * t40 - t54 * t39) * MDP(17) + (g(1) * t39 - t54 * t40) * MDP(18);];
taug = t1;
