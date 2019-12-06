% Calculate Gravitation load on the joints for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:59
% EndTime: 2019-12-05 17:07:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (139->22), mult. (95->29), div. (0->0), fcn. (74->8), ass. (0->15)
t36 = qJ(4) + qJ(5);
t33 = sin(t36);
t34 = cos(t36);
t35 = pkin(9) + qJ(2);
t32 = qJ(3) + t35;
t28 = sin(t32);
t29 = cos(t32);
t41 = g(1) * t29 + g(2) * t28;
t42 = (-g(3) * t34 + t41 * t33) * MDP(20) + (g(3) * t33 + t41 * t34) * MDP(21);
t37 = sin(qJ(4));
t38 = cos(qJ(4));
t39 = t41 * MDP(7) + (t38 * MDP(13) - t37 * MDP(14) + t34 * MDP(20) - t33 * MDP(21) + MDP(6)) * (g(1) * t28 - g(2) * t29);
t31 = cos(t35);
t30 = sin(t35);
t1 = [-g(3) * MDP(1); (g(1) * t30 - g(2) * t31) * MDP(3) + (g(1) * t31 + g(2) * t30) * MDP(4) + t39; t39; (-g(3) * t38 + t41 * t37) * MDP(13) + (g(3) * t37 + t41 * t38) * MDP(14) + t42; t42;];
taug = t1;
