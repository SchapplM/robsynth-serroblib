% Calculate Gravitation load on the joints for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (105->23), mult. (99->36), div. (0->0), fcn. (88->8), ass. (0->18)
t31 = pkin(9) + qJ(3);
t30 = qJ(4) + t31;
t26 = sin(t30);
t44 = g(3) * t26;
t32 = sin(pkin(8));
t34 = sin(qJ(5));
t42 = t32 * t34;
t35 = cos(qJ(5));
t41 = t32 * t35;
t33 = cos(pkin(8));
t40 = t33 * t34;
t39 = t33 * t35;
t27 = cos(t30);
t37 = g(1) * t33 + g(2) * t32;
t38 = (t37 * t27 + t44) * MDP(8) + (t35 * MDP(14) - t34 * MDP(15) + MDP(7)) * (-g(3) * t27 + t37 * t26);
t29 = cos(t31);
t28 = sin(t31);
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t32 + g(2) * t33) * MDP(2); (-g(3) * t29 + t37 * t28) * MDP(4) + (g(3) * t28 + t37 * t29) * MDP(5) + t38; t38; (-g(1) * (-t27 * t40 + t41) - g(2) * (-t27 * t42 - t39) + t34 * t44) * MDP(14) + (-g(1) * (-t27 * t39 - t42) - g(2) * (-t27 * t41 + t40) + t35 * t44) * MDP(15);];
taug = t1;
