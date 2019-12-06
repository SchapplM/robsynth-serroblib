% Calculate Gravitation load on the joints for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:20
% EndTime: 2019-12-05 15:01:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (94->28), mult. (105->40), div. (0->0), fcn. (90->8), ass. (0->17)
t33 = sin(pkin(7));
t35 = cos(pkin(7));
t38 = g(1) * t35 + g(2) * t33;
t31 = pkin(8) + qJ(3);
t27 = sin(t31);
t29 = cos(t31);
t23 = -g(3) * t29 + t38 * t27;
t45 = g(3) * t27;
t30 = pkin(9) + qJ(5);
t26 = sin(t30);
t43 = t33 * t26;
t28 = cos(t30);
t42 = t33 * t28;
t41 = t35 * t26;
t40 = t35 * t28;
t39 = MDP(2) + MDP(9);
t1 = [(-MDP(1) - t39) * g(3); t39 * (-g(1) * t33 + g(2) * t35); (-g(3) * (t29 * pkin(3) + t27 * qJ(4)) + t38 * (pkin(3) * t27 - qJ(4) * t29)) * MDP(9) + (MDP(5) - MDP(8)) * (t38 * t29 + t45) + (MDP(15) * t28 - MDP(16) * t26 + MDP(6) * cos(pkin(9)) - MDP(7) * sin(pkin(9)) + MDP(4)) * t23; -t23 * MDP(9); (-g(1) * (-t29 * t41 + t42) - g(2) * (-t29 * t43 - t40) + t26 * t45) * MDP(15) + (-g(1) * (-t29 * t40 - t43) - g(2) * (-t29 * t42 + t41) + t28 * t45) * MDP(16);];
taug = t1;
