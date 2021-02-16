% Calculate Gravitation load on the joints for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:17
% EndTime: 2021-01-15 14:48:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (108->30), mult. (146->40), div. (0->0), fcn. (129->6), ass. (0->20)
t55 = MDP(11) + MDP(13);
t54 = MDP(12) + MDP(14);
t38 = sin(pkin(7));
t39 = cos(pkin(7));
t44 = g(1) * t39 + g(2) * t38;
t37 = pkin(8) + qJ(3);
t35 = sin(t37);
t36 = cos(t37);
t31 = -g(3) * t36 + t35 * t44;
t51 = g(3) * t35;
t41 = sin(qJ(4));
t49 = t38 * t41;
t42 = cos(qJ(4));
t48 = t38 * t42;
t47 = t39 * t41;
t46 = t39 * t42;
t45 = MDP(16) + MDP(2);
t40 = -qJ(5) - pkin(6);
t34 = pkin(4) * t42 + pkin(3);
t1 = [(-MDP(1) - t45) * g(3); t45 * (-g(1) * t38 + g(2) * t39); (-g(3) * (t34 * t36 - t35 * t40) + t44 * (t34 * t35 + t36 * t40)) * MDP(16) + (MDP(5) - MDP(15)) * (t36 * t44 + t51) + (-t41 * t54 + t42 * t55 + MDP(4)) * t31; t54 * (-g(1) * (-t36 * t46 - t49) - g(2) * (-t36 * t48 + t47) + t42 * t51) + (MDP(16) * pkin(4) + t55) * (-g(1) * (-t36 * t47 + t48) - g(2) * (-t36 * t49 - t46) + t41 * t51); -t31 * MDP(16);];
taug = t1;
