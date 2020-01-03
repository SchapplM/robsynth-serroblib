% Calculate Gravitation load on the joints for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (102->29), mult. (120->49), div. (0->0), fcn. (101->10), ass. (0->21)
t36 = qJ(2) + qJ(3);
t33 = pkin(9) + t36;
t31 = sin(t33);
t51 = g(3) * t31;
t37 = sin(pkin(8));
t39 = sin(qJ(5));
t49 = t37 * t39;
t41 = cos(qJ(5));
t48 = t37 * t41;
t38 = cos(pkin(8));
t47 = t38 * t39;
t46 = t38 * t41;
t32 = cos(t33);
t34 = sin(t36);
t35 = cos(t36);
t44 = g(1) * t38 + g(2) * t37;
t43 = -g(3) * t35 + t44 * t34;
t45 = t43 * MDP(6) + (g(3) * t34 + t44 * t35) * MDP(7) + (t41 * MDP(14) - t39 * MDP(15)) * (-g(3) * t32 + t44 * t31);
t42 = cos(qJ(2));
t40 = sin(qJ(2));
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(3) * t42 + t44 * t40) * MDP(3) + (g(3) * t40 + t44 * t42) * MDP(4) + (-g(3) * (t42 * pkin(2) + pkin(3) * t35) - t44 * (-t40 * pkin(2) - pkin(3) * t34)) * MDP(8) + t45; t43 * MDP(8) * pkin(3) + t45; (-g(1) * t37 + g(2) * t38) * MDP(8); (-g(1) * (-t32 * t47 + t48) - g(2) * (-t32 * t49 - t46) + t39 * t51) * MDP(14) + (-g(1) * (-t32 * t46 - t49) - g(2) * (-t32 * t48 + t47) + t41 * t51) * MDP(15);];
taug = t1;
