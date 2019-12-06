% Calculate Gravitation load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:22
% EndTime: 2019-12-05 15:03:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->26), mult. (93->38), div. (0->0), fcn. (79->6), ass. (0->16)
t30 = sin(pkin(7));
t31 = cos(pkin(7));
t36 = g(1) * t31 + g(2) * t30;
t29 = pkin(8) + qJ(3);
t28 = cos(t29);
t42 = g(3) * t28;
t32 = sin(qJ(5));
t41 = t30 * t32;
t33 = cos(qJ(5));
t40 = t30 * t33;
t39 = t31 * t32;
t38 = t31 * t33;
t37 = MDP(2) + MDP(8);
t27 = sin(t29);
t23 = t36 * t27 - t42;
t1 = [(-MDP(1) - t37) * g(3); t37 * (-g(1) * t30 + g(2) * t31); (-g(3) * (t28 * pkin(3) + t27 * qJ(4)) + t36 * (pkin(3) * t27 - qJ(4) * t28)) * MDP(8) + (MDP(4) - MDP(6)) * t23 + (-MDP(14) * t32 - MDP(15) * t33 + MDP(5) - MDP(7)) * (g(3) * t27 + t36 * t28); -t23 * MDP(8); (-g(1) * (t27 * t38 - t41) - g(2) * (t27 * t40 + t39) + t33 * t42) * MDP(14) + (-g(1) * (-t27 * t39 - t40) - g(2) * (-t27 * t41 + t38) - t32 * t42) * MDP(15);];
taug = t1;
