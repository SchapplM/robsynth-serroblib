% Calculate Gravitation load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:35
% EndTime: 2019-12-05 18:14:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (141->20), mult. (84->28), div. (0->0), fcn. (60->8), ass. (0->14)
t29 = qJ(1) + pkin(9) + qJ(3);
t28 = qJ(4) + t29;
t24 = sin(t28);
t25 = cos(t28);
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t36 = -g(2) * t24 + g(3) * t25;
t38 = t36 * MDP(10) + (t32 * MDP(16) - t30 * MDP(17) + MDP(9)) * (g(2) * t25 + g(3) * t24);
t26 = sin(t29);
t27 = cos(t29);
t34 = (g(2) * t27 + g(3) * t26) * MDP(6) + (-g(2) * t26 + g(3) * t27) * MDP(7) + t38;
t33 = cos(qJ(1));
t31 = sin(qJ(1));
t1 = [(-g(2) * t31 + g(3) * t33) * MDP(3) + t34 + (MDP(4) * pkin(1) + MDP(2)) * (g(2) * t33 + g(3) * t31); -g(1) * MDP(4); t34; t38; (-g(1) * t32 + t36 * t30) * MDP(16) + (g(1) * t30 + t36 * t32) * MDP(17);];
taug = t1;
