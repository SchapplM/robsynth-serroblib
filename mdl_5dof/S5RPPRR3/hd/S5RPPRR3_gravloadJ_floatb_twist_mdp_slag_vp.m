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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:15
% EndTime: 2020-01-03 11:28:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (117->30), mult. (94->42), div. (0->0), fcn. (70->10), ass. (0->15)
t28 = pkin(9) + qJ(4);
t27 = qJ(5) + t28;
t21 = sin(t27);
t22 = cos(t27);
t29 = qJ(1) + pkin(8);
t24 = sin(t29);
t26 = cos(t29);
t35 = g(2) * t24 - g(3) * t26;
t37 = (-g(1) * t22 + t35 * t21) * MDP(21) + (g(1) * t21 + t35 * t22) * MDP(22);
t36 = g(2) * t26 + g(3) * t24;
t33 = cos(qJ(1));
t32 = sin(qJ(1));
t25 = cos(t28);
t23 = sin(t28);
t1 = [(g(2) * t32 - g(3) * t33) * MDP(3) - t35 * MDP(7) + (-g(2) * (t33 * pkin(1) + t26 * pkin(2) + t24 * qJ(3)) - g(3) * (t32 * pkin(1) + t24 * pkin(2) - t26 * qJ(3))) * MDP(8) + (pkin(1) * MDP(4) + MDP(2)) * (-g(2) * t33 - g(3) * t32) + (-t25 * MDP(14) + t23 * MDP(15) - MDP(21) * t22 + MDP(22) * t21 - MDP(5) * cos(pkin(9)) + MDP(6) * sin(pkin(9))) * t36; (-MDP(4) - MDP(8)) * g(1); t36 * MDP(8); (-g(1) * t25 + t35 * t23) * MDP(14) + (g(1) * t23 + t35 * t25) * MDP(15) + t37; t37;];
taug = t1;
