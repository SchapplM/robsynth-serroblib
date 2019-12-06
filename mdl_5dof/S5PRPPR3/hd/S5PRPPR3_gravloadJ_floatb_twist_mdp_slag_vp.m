% Calculate Gravitation load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (72->30), mult. (107->47), div. (0->0), fcn. (89->8), ass. (0->18)
t28 = sin(pkin(7));
t29 = cos(pkin(7));
t37 = g(1) * t29 + g(2) * t28;
t27 = qJ(2) + pkin(8);
t26 = cos(t27);
t43 = g(3) * t26;
t30 = sin(qJ(5));
t42 = t28 * t30;
t32 = cos(qJ(5));
t41 = t28 * t32;
t40 = t29 * t30;
t39 = t29 * t32;
t38 = MDP(5) + MDP(8);
t33 = cos(qJ(2));
t31 = sin(qJ(2));
t25 = sin(t27);
t23 = -t37 * t25 + t43;
t1 = [(-MDP(1) - t38) * g(3); (g(3) * t31 + t37 * t33) * MDP(4) + t23 * MDP(6) + (-g(3) * (t33 * pkin(2) + t26 * pkin(3) + t25 * qJ(4)) + t37 * (pkin(2) * t31 + pkin(3) * t25 - qJ(4) * t26)) * MDP(8) + (pkin(2) * MDP(5) + MDP(3)) * (-g(3) * t33 + t37 * t31) + (MDP(14) * t30 + MDP(15) * t32 + MDP(7)) * (-g(3) * t25 - t37 * t26); t38 * (-g(1) * t28 + g(2) * t29); t23 * MDP(8); (-g(1) * (t25 * t39 - t42) - g(2) * (t25 * t41 + t40) + t32 * t43) * MDP(14) + (-g(1) * (-t25 * t40 - t41) - g(2) * (-t25 * t42 + t39) - t30 * t43) * MDP(15);];
taug = t1;
