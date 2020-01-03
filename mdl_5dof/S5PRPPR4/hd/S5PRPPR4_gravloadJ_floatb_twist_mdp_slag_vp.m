% Calculate Gravitation load on the joints for
% S5PRPPR4
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:56
% EndTime: 2019-12-31 17:36:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (103->34), mult. (119->49), div. (0->0), fcn. (105->6), ass. (0->21)
t39 = pkin(7) + qJ(2);
t37 = sin(t39);
t49 = g(1) * t37;
t38 = cos(t39);
t48 = t38 * pkin(2) + t37 * qJ(3);
t47 = -MDP(12) - MDP(8);
t32 = g(1) * t38 + g(2) * t37;
t31 = -g(2) * t38 + t49;
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t46 = pkin(3) * t41 + qJ(4) * t40;
t42 = sin(qJ(5));
t43 = cos(qJ(5));
t45 = t40 * t43 - t41 * t42;
t44 = t40 * t42 + t41 * t43;
t34 = t38 * qJ(3);
t28 = t44 * t38;
t27 = t45 * t38;
t26 = t44 * t37;
t25 = t45 * t37;
t1 = [(-MDP(1) + t47) * g(3); (-g(1) * (-t37 * pkin(2) + t34) - g(2) * t48) * MDP(8) + (-g(1) * t34 - g(2) * (t46 * t38 + t48) - (-pkin(2) - t46) * t49) * MDP(12) + (g(1) * t26 - g(2) * t28) * MDP(18) + (g(1) * t25 - g(2) * t27) * MDP(19) + (MDP(4) - MDP(7) - MDP(10)) * t32 + (MDP(3) + (MDP(5) + MDP(9)) * t41 + (-MDP(6) + MDP(11)) * t40) * t31; t47 * t31; (g(3) * t41 - t32 * t40) * MDP(12); (-g(1) * t27 - g(2) * t25 + g(3) * t44) * MDP(18) + (g(1) * t28 + g(2) * t26 + g(3) * t45) * MDP(19);];
taug = t1;
