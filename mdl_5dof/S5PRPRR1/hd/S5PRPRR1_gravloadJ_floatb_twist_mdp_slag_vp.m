% Calculate Gravitation load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:05
% DurationCPUTime: 0.07s
% Computational Cost: add. (117->26), mult. (87->33), div. (0->0), fcn. (66->8), ass. (0->13)
t27 = pkin(9) + qJ(4);
t26 = qJ(5) + t27;
t20 = sin(t26);
t21 = cos(t26);
t28 = pkin(8) + qJ(2);
t23 = sin(t28);
t25 = cos(t28);
t31 = g(1) * t25 + g(2) * t23;
t32 = (-g(3) * t21 + t31 * t20) * MDP(21) + (g(3) * t20 + t31 * t21) * MDP(22);
t18 = g(1) * t23 - g(2) * t25;
t24 = cos(t27);
t22 = sin(t27);
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(1) * (-t23 * pkin(2) + t25 * qJ(3)) - g(2) * (t25 * pkin(2) + t23 * qJ(3))) * MDP(8) + (MDP(4) - MDP(7)) * t31 + (t24 * MDP(14) - t22 * MDP(15) + MDP(21) * t21 - MDP(22) * t20 + MDP(5) * cos(pkin(9)) - MDP(6) * sin(pkin(9)) + MDP(3)) * t18; -t18 * MDP(8); (-g(3) * t24 + t31 * t22) * MDP(14) + (g(3) * t22 + t31 * t24) * MDP(15) + t32; t32;];
taug = t1;
