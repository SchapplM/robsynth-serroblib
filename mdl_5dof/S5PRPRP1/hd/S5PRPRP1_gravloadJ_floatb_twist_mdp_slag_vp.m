% Calculate Gravitation load on the joints for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:40
% EndTime: 2019-12-05 15:28:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (142->34), mult. (119->41), div. (0->0), fcn. (89->6), ass. (0->17)
t36 = pkin(7) + qJ(2);
t32 = sin(t36);
t34 = cos(t36);
t28 = g(1) * t34 + g(2) * t32;
t47 = MDP(14) + MDP(16);
t46 = MDP(15) - MDP(18);
t43 = -MDP(19) - MDP(8);
t27 = g(1) * t32 - g(2) * t34;
t35 = pkin(8) + qJ(4);
t31 = sin(t35);
t33 = cos(t35);
t42 = pkin(4) * t33 + qJ(5) * t31;
t38 = cos(pkin(8));
t40 = pkin(3) * t38 + pkin(2) + t42;
t39 = -pkin(6) - qJ(3);
t23 = -g(3) * t33 + t28 * t31;
t1 = [(-MDP(1) + t43) * g(3); (-g(1) * (-pkin(2) * t32 + qJ(3) * t34) - g(2) * (pkin(2) * t34 + qJ(3) * t32)) * MDP(8) + ((g(1) * t39 - g(2) * t40) * t34 + (g(1) * t40 + g(2) * t39) * t32) * MDP(19) + (MDP(4) - MDP(7) - MDP(17)) * t28 + (MDP(5) * t38 - MDP(6) * sin(pkin(8)) - t46 * t31 + t47 * t33 + MDP(3)) * t27; t43 * t27; (-g(3) * t42 + t28 * (pkin(4) * t31 - qJ(5) * t33)) * MDP(19) + t46 * (g(3) * t31 + t28 * t33) + t47 * t23; -t23 * MDP(19);];
taug = t1;
