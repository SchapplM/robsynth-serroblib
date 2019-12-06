% Calculate Gravitation load on the joints for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (102->23), mult. (108->37), div. (0->0), fcn. (93->8), ass. (0->17)
t28 = qJ(2) + pkin(9) + qJ(4);
t26 = sin(t28);
t44 = g(3) * t26;
t29 = sin(pkin(8));
t31 = sin(qJ(5));
t42 = t29 * t31;
t33 = cos(qJ(5));
t41 = t29 * t33;
t30 = cos(pkin(8));
t40 = t30 * t31;
t39 = t30 * t33;
t27 = cos(t28);
t37 = g(1) * t30 + g(2) * t29;
t38 = (t37 * t27 + t44) * MDP(8) + (t33 * MDP(14) - t31 * MDP(15) + MDP(7)) * (-g(3) * t27 + t37 * t26);
t34 = cos(qJ(2));
t32 = sin(qJ(2));
t1 = [(-MDP(1) - MDP(5)) * g(3); (g(3) * t32 + t37 * t34) * MDP(4) + t38 + (MDP(5) * pkin(2) + MDP(3)) * (-g(3) * t34 + t37 * t32); (-g(1) * t29 + g(2) * t30) * MDP(5); t38; (-g(1) * (-t27 * t40 + t41) - g(2) * (-t27 * t42 - t39) + t31 * t44) * MDP(14) + (-g(1) * (-t27 * t39 - t42) - g(2) * (-t27 * t41 + t40) + t33 * t44) * MDP(15);];
taug = t1;
