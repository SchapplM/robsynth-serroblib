% Calculate Gravitation load on the joints for
% S5RRPPR2
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
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:22
% EndTime: 2019-12-05 18:20:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (177->44), mult. (139->66), div. (0->0), fcn. (118->10), ass. (0->29)
t52 = qJ(1) + qJ(2);
t50 = sin(t52);
t70 = pkin(2) * t50;
t51 = cos(t52);
t69 = pkin(2) * t51;
t53 = sin(pkin(9));
t68 = g(1) * t53;
t56 = sin(qJ(1));
t67 = t56 * pkin(1);
t58 = cos(qJ(1));
t66 = t58 * pkin(1);
t54 = cos(pkin(9));
t55 = sin(qJ(5));
t65 = t54 * t55;
t57 = cos(qJ(5));
t64 = t54 * t57;
t49 = pkin(8) + t52;
t47 = sin(t49);
t48 = cos(t49);
t63 = g(2) * t48 + g(3) * t47;
t62 = g(2) * t51 + g(3) * t50;
t61 = -t47 * pkin(3) + t48 * qJ(4) - t70;
t39 = t47 * t65 + t48 * t57;
t40 = t47 * t64 - t48 * t55;
t41 = -t47 * t57 + t48 * t65;
t42 = -t47 * t55 - t48 * t64;
t60 = (-g(2) * t41 - g(3) * t39) * MDP(18) + (-g(2) * t42 + g(3) * t40) * MDP(17) + (g(2) * t47 - g(3) * t48) * MDP(10) + (-g(2) * t50 + g(3) * t51) * MDP(6) + t62 * MDP(5) + (t54 * MDP(8) - t53 * MDP(9)) * t63;
t59 = -t48 * pkin(3) - t47 * qJ(4) - t69;
t1 = [(g(2) * t58 + g(3) * t56) * MDP(2) + (-g(2) * t56 + g(3) * t58) * MDP(3) + (-g(2) * (-t66 - t69) - g(3) * (-t67 - t70)) * MDP(7) + (-g(2) * (t59 - t66) - g(3) * (t61 - t67)) * MDP(11) + t60; t62 * pkin(2) * MDP(7) + (-g(2) * t59 - g(3) * t61) * MDP(11) + t60; (-MDP(11) - MDP(7)) * g(1); -t63 * MDP(11); (-g(2) * t39 + g(3) * t41 + t55 * t68) * MDP(17) + (-g(2) * t40 - g(3) * t42 + t57 * t68) * MDP(18);];
taug = t1;
