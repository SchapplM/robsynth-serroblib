% Calculate Gravitation load on the joints for
% S5RRRPR1
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:51
% EndTime: 2019-12-05 18:38:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (162->35), mult. (146->47), div. (0->0), fcn. (113->8), ass. (0->19)
t39 = qJ(2) + qJ(3);
t34 = pkin(9) + qJ(5) + t39;
t31 = sin(t34);
t32 = cos(t34);
t41 = sin(qJ(1));
t43 = cos(qJ(1));
t45 = g(1) * t43 + g(2) * t41;
t48 = (-g(3) * t32 + t45 * t31) * MDP(25) + (g(3) * t31 + t45 * t32) * MDP(26);
t36 = cos(t39);
t42 = cos(qJ(2));
t47 = t42 * pkin(2) + pkin(3) * t36;
t35 = sin(t39);
t44 = -g(3) * t36 + t45 * t35;
t46 = t44 * MDP(16) + (g(3) * t35 + t45 * t36) * MDP(17) + t48;
t29 = g(1) * t41 - g(2) * t43;
t40 = sin(qJ(2));
t38 = -qJ(4) - pkin(7) - pkin(6);
t27 = pkin(1) + t47;
t1 = [(-g(1) * (-t41 * t27 - t43 * t38) - g(2) * (t43 * t27 - t41 * t38)) * MDP(19) + (MDP(3) - MDP(18)) * t45 + (-t40 * MDP(10) + MDP(16) * t36 - MDP(17) * t35 + MDP(25) * t32 - MDP(26) * t31 + t42 * MDP(9) + MDP(2)) * t29; (-g(3) * t42 + t45 * t40) * MDP(9) + (g(3) * t40 + t45 * t42) * MDP(10) + (-g(3) * t47 - t45 * (-t40 * pkin(2) - pkin(3) * t35)) * MDP(19) + t46; t44 * MDP(19) * pkin(3) + t46; -t29 * MDP(19); t48;];
taug = t1;
