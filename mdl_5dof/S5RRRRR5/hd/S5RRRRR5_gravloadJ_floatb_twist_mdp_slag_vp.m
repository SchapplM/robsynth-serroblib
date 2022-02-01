% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:07
% EndTime: 2022-01-20 12:02:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (189->25), mult. (132->34), div. (0->0), fcn. (102->10), ass. (0->18)
t41 = qJ(4) + qJ(5);
t36 = sin(t41);
t38 = cos(t41);
t42 = qJ(1) + qJ(2);
t40 = qJ(3) + t42;
t34 = sin(t40);
t35 = cos(t40);
t50 = g(1) * t35 + g(2) * t34;
t51 = (-g(3) * t38 + t50 * t36) * MDP(22) + (g(3) * t36 + t50 * t38) * MDP(23);
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t48 = t50 * MDP(9) + (t45 * MDP(15) - t43 * MDP(16) + t38 * MDP(22) - t36 * MDP(23) + MDP(8)) * (g(1) * t34 - g(2) * t35);
t37 = sin(t42);
t39 = cos(t42);
t47 = (g(1) * t37 - g(2) * t39) * MDP(5) + (g(1) * t39 + g(2) * t37) * MDP(6) + t48;
t46 = cos(qJ(1));
t44 = sin(qJ(1));
t1 = [(g(1) * t44 - g(2) * t46) * MDP(2) + (g(1) * t46 + g(2) * t44) * MDP(3) + t47; t47; t48; (-g(3) * t45 + t50 * t43) * MDP(15) + (g(3) * t43 + t50 * t45) * MDP(16) + t51; t51;];
taug = t1;
