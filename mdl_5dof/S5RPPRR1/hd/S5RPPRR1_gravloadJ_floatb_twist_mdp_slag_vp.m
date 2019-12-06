% Calculate Gravitation load on the joints for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (65->30), mult. (99->36), div. (0->0), fcn. (74->6), ass. (0->13)
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t20 = g(1) * t31 + g(2) * t29;
t27 = qJ(4) + qJ(5);
t21 = sin(t27);
t22 = cos(t27);
t33 = (g(3) * t21 - t20 * t22) * MDP(22) + (g(3) * t22 + t20 * t21) * MDP(23);
t32 = t31 * pkin(1) + t29 * qJ(2);
t19 = g(1) * t29 - g(2) * t31;
t30 = cos(qJ(4));
t28 = sin(qJ(4));
t24 = t31 * qJ(2);
t1 = [(-g(1) * (-t29 * pkin(1) + t24) - g(2) * t32) * MDP(6) + (-g(1) * (t24 + (-pkin(1) - qJ(3)) * t29) - g(2) * (t31 * qJ(3) + t32)) * MDP(9) + (MDP(3) - MDP(5) - MDP(7)) * t20 + (t28 * MDP(15) + t30 * MDP(16) + MDP(22) * t21 + MDP(23) * t22 + MDP(2) - MDP(4) + MDP(8)) * t19; (-MDP(6) - MDP(9)) * t19; -t20 * MDP(9); (g(3) * t28 - t20 * t30) * MDP(15) + (g(3) * t30 + t20 * t28) * MDP(16) + t33; t33;];
taug = t1;
