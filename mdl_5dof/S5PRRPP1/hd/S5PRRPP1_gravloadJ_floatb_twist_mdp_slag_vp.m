% Calculate Gravitation load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:54
% EndTime: 2019-12-05 16:06:55
% DurationCPUTime: 0.16s
% Computational Cost: add. (138->40), mult. (125->50), div. (0->0), fcn. (93->6), ass. (0->18)
t33 = pkin(7) + qJ(2);
t28 = sin(t33);
t30 = cos(t33);
t24 = g(1) * t30 + g(2) * t28;
t41 = -MDP(13) - MDP(17);
t23 = g(1) * t28 - g(2) * t30;
t34 = qJ(3) + pkin(8);
t29 = sin(t34);
t31 = cos(t34);
t40 = t31 * pkin(4) + t29 * qJ(5);
t37 = cos(qJ(3));
t36 = sin(qJ(3));
t35 = -qJ(4) - pkin(6);
t32 = t37 * pkin(3);
t27 = t32 + pkin(2);
t25 = t30 * t27;
t22 = -g(3) * t31 + t24 * t29;
t1 = [(-MDP(1) + t41) * g(3); (-g(1) * (-t28 * t27 - t30 * t35) - g(2) * (-t28 * t35 + t25)) * MDP(13) + (-g(2) * t25 + (g(1) * t35 - g(2) * t40) * t30 + (-g(1) * (-t27 - t40) + g(2) * t35) * t28) * MDP(17) + (MDP(4) - MDP(12) - MDP(15)) * t24 + (t37 * MDP(10) - t36 * MDP(11) + t31 * MDP(14) + t29 * MDP(16) + MDP(3)) * t23; (g(3) * t36 + t24 * t37) * MDP(11) + t22 * MDP(14) + (-g(3) * t29 - t24 * t31) * MDP(16) + (-g(3) * (t32 + t40) + t24 * (pkin(3) * t36 + pkin(4) * t29 - qJ(5) * t31)) * MDP(17) + (pkin(3) * MDP(13) + MDP(10)) * (-g(3) * t37 + t24 * t36); t41 * t23; -t22 * MDP(17);];
taug = t1;
