% Calculate Gravitation load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (113->30), mult. (89->41), div. (0->0), fcn. (66->9), ass. (0->15)
t27 = pkin(9) + qJ(4);
t26 = qJ(5) + t27;
t20 = sin(t26);
t21 = cos(t26);
t28 = qJ(1) + pkin(8);
t23 = sin(t28);
t25 = cos(t28);
t34 = g(1) * t25 + g(2) * t23;
t35 = (-g(3) * t21 + t34 * t20) * MDP(20) + (g(3) * t20 + t34 * t21) * MDP(21);
t33 = g(1) * t23 - g(2) * t25;
t31 = cos(qJ(1));
t30 = sin(qJ(1));
t24 = cos(t27);
t22 = sin(t27);
t1 = [(g(1) * t31 + g(2) * t30) * MDP(3) - t34 * MDP(6) + (-g(1) * (-t30 * pkin(1) - t23 * pkin(2) + t25 * qJ(3)) - g(2) * (t31 * pkin(1) + t25 * pkin(2) + t23 * qJ(3))) * MDP(7) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t30 - g(2) * t31) + (t24 * MDP(13) - t22 * MDP(14) + MDP(20) * t21 - MDP(21) * t20 + MDP(5) * cos(pkin(9))) * t33; (-MDP(4) - MDP(7)) * g(3); -t33 * MDP(7); (-g(3) * t24 + t34 * t22) * MDP(13) + (g(3) * t22 + t34 * t24) * MDP(14) + t35; t35;];
taug = t1;
