% Calculate Gravitation load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:59
% EndTime: 2019-12-05 18:40:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (147->31), mult. (100->46), div. (0->0), fcn. (70->10), ass. (0->20)
t38 = qJ(1) + qJ(2);
t37 = qJ(3) + t38;
t32 = pkin(9) + t37;
t30 = sin(t32);
t31 = cos(t32);
t33 = sin(t37);
t34 = cos(t37);
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t44 = g(2) * t34 + g(3) * t33;
t49 = (-g(2) * t33 + g(3) * t34) * MDP(9) + t44 * MDP(8) + (t41 * MDP(16) - t39 * MDP(17)) * (g(2) * t31 + g(3) * t30);
t35 = sin(t38);
t48 = -pkin(2) * t35 - pkin(3) * t33;
t36 = cos(t38);
t47 = -pkin(2) * t36 - pkin(3) * t34;
t45 = -g(2) * t30 + g(3) * t31;
t43 = (-g(2) * t35 + g(3) * t36) * MDP(6) + (g(2) * t36 + g(3) * t35) * MDP(5) + t49;
t42 = cos(qJ(1));
t40 = sin(qJ(1));
t1 = [(g(2) * t42 + g(3) * t40) * MDP(2) + (-g(2) * t40 + g(3) * t42) * MDP(3) + (-g(2) * (-t42 * pkin(1) + t47) - g(3) * (-t40 * pkin(1) + t48)) * MDP(10) + t43; (-g(2) * t47 - g(3) * t48) * MDP(10) + t43; t44 * MDP(10) * pkin(3) + t49; -g(1) * MDP(10); (-g(1) * t41 + t45 * t39) * MDP(16) + (g(1) * t39 + t45 * t41) * MDP(17);];
taug = t1;
