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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (147->31), mult. (100->46), div. (0->0), fcn. (70->10), ass. (0->20)
t42 = qJ(1) + qJ(2);
t41 = qJ(3) + t42;
t38 = cos(t41);
t40 = cos(t42);
t53 = pkin(2) * t40 + pkin(3) * t38;
t36 = pkin(9) + t41;
t32 = sin(t36);
t33 = cos(t36);
t37 = sin(t41);
t43 = sin(qJ(5));
t45 = cos(qJ(5));
t48 = g(1) * t37 - g(2) * t38;
t52 = t48 * MDP(8) + (g(1) * t38 + g(2) * t37) * MDP(9) + (t45 * MDP(16) - t43 * MDP(17)) * (g(1) * t32 - g(2) * t33);
t39 = sin(t42);
t51 = -pkin(2) * t39 - pkin(3) * t37;
t50 = g(1) * t33 + g(2) * t32;
t47 = (g(1) * t39 - g(2) * t40) * MDP(5) + (g(1) * t40 + g(2) * t39) * MDP(6) + t52;
t46 = cos(qJ(1));
t44 = sin(qJ(1));
t1 = [(g(1) * t44 - g(2) * t46) * MDP(2) + (g(1) * t46 + g(2) * t44) * MDP(3) + (-g(1) * (-t44 * pkin(1) + t51) - g(2) * (t46 * pkin(1) + t53)) * MDP(10) + t47; (-g(1) * t51 - g(2) * t53) * MDP(10) + t47; t48 * MDP(10) * pkin(3) + t52; -g(3) * MDP(10); (-g(3) * t45 + t50 * t43) * MDP(16) + (g(3) * t43 + t50 * t45) * MDP(17);];
taug = t1;
