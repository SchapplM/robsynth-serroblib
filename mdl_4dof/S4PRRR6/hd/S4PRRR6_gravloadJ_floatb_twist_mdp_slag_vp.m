% Calculate Gravitation load on the joints for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:07
% EndTime: 2019-12-31 16:35:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->29), mult. (115->50), div. (0->0), fcn. (114->8), ass. (0->17)
t25 = qJ(3) + qJ(4);
t23 = sin(t25);
t24 = cos(t25);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t43 = t30 * MDP(10) - t28 * MDP(11) + t24 * MDP(17) - t23 * MDP(18) + MDP(3);
t29 = sin(qJ(2));
t42 = g(3) * t29;
t26 = sin(pkin(7));
t31 = cos(qJ(2));
t41 = t26 * t31;
t27 = cos(pkin(7));
t40 = t27 * t31;
t39 = t28 * t31;
t38 = t30 * t31;
t37 = (-g(1) * (-t23 * t40 + t26 * t24) - g(2) * (-t23 * t41 - t27 * t24) + t23 * t42) * MDP(17) + (-g(1) * (-t26 * t23 - t24 * t40) - g(2) * (t27 * t23 - t24 * t41) + t24 * t42) * MDP(18);
t1 = [-g(3) * MDP(1); (t29 * MDP(4) - t43 * t31) * g(3) + (MDP(4) * t31 + t43 * t29) * (g(1) * t27 + g(2) * t26); (-g(1) * (t26 * t30 - t27 * t39) - g(2) * (-t26 * t39 - t27 * t30) + t28 * t42) * MDP(10) + (-g(1) * (-t26 * t28 - t27 * t38) - g(2) * (-t26 * t38 + t27 * t28) + t30 * t42) * MDP(11) + t37; t37;];
taug = t1;
