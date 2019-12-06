% Calculate Gravitation load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:19
% EndTime: 2019-12-05 17:06:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (143->20), mult. (79->27), div. (0->0), fcn. (58->8), ass. (0->15)
t32 = pkin(9) + qJ(2);
t31 = qJ(3) + t32;
t28 = qJ(4) + t31;
t24 = sin(t28);
t25 = cos(t28);
t33 = sin(qJ(5));
t34 = cos(qJ(5));
t37 = g(1) * t25 + g(2) * t24;
t38 = t37 * MDP(10) + (t34 * MDP(16) - t33 * MDP(17) + MDP(9)) * (g(1) * t24 - g(2) * t25);
t26 = sin(t31);
t27 = cos(t31);
t35 = (g(1) * t26 - g(2) * t27) * MDP(6) + (g(1) * t27 + g(2) * t26) * MDP(7) + t38;
t30 = cos(t32);
t29 = sin(t32);
t1 = [-g(3) * MDP(1); (g(1) * t29 - g(2) * t30) * MDP(3) + (g(1) * t30 + g(2) * t29) * MDP(4) + t35; t35; t38; (-g(3) * t34 + t37 * t33) * MDP(16) + (g(3) * t33 + t37 * t34) * MDP(17);];
taug = t1;
