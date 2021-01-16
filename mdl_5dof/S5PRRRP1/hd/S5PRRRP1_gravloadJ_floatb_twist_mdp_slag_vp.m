% Calculate Gravitation load on the joints for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:41
% EndTime: 2021-01-15 16:14:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (166->30), mult. (122->36), div. (0->0), fcn. (89->6), ass. (0->18)
t51 = MDP(13) + MDP(15);
t50 = MDP(14) + MDP(16);
t42 = pkin(8) + qJ(2);
t41 = qJ(3) + t42;
t36 = sin(t41);
t37 = cos(t41);
t45 = cos(qJ(4));
t38 = t45 * pkin(4) + pkin(3);
t43 = qJ(5) + pkin(7);
t49 = t36 * t43 + t37 * t38;
t48 = -t36 * t38 + t43 * t37;
t33 = g(1) * t37 + g(2) * t36;
t32 = g(1) * t36 - g(2) * t37;
t44 = sin(qJ(4));
t47 = (-MDP(17) + MDP(7)) * t33 + (-t50 * t44 + t51 * t45 + MDP(6)) * t32;
t40 = cos(t42);
t39 = sin(t42);
t1 = [(-MDP(1) - MDP(18)) * g(3); (g(1) * t39 - g(2) * t40) * MDP(3) + (g(1) * t40 + g(2) * t39) * MDP(4) + (-g(1) * (-pkin(2) * t39 + t48) - g(2) * (pkin(2) * t40 + t49)) * MDP(18) + t47; (-g(1) * t48 - g(2) * t49) * MDP(18) + t47; t50 * (g(3) * t44 + t33 * t45) + (MDP(18) * pkin(4) + t51) * (-g(3) * t45 + t33 * t44); -t32 * MDP(18);];
taug = t1;
