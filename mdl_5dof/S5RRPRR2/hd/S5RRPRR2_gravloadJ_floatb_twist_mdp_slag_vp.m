% Calculate Gravitation load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:27
% DurationCPUTime: 0.08s
% Computational Cost: add. (160->31), mult. (132->40), div. (0->0), fcn. (103->8), ass. (0->17)
t32 = qJ(2) + pkin(9) + qJ(4);
t30 = qJ(5) + t32;
t26 = sin(t30);
t27 = cos(t30);
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t39 = g(1) * t37 + g(2) * t35;
t41 = (-g(3) * t27 + t39 * t26) * MDP(25) + (g(3) * t26 + t39 * t27) * MDP(26);
t28 = sin(t32);
t29 = cos(t32);
t40 = (-g(3) * t29 + t39 * t28) * MDP(18) + (g(3) * t28 + t39 * t29) * MDP(19) + t41;
t24 = g(1) * t35 - g(2) * t37;
t36 = cos(qJ(2));
t34 = sin(qJ(2));
t33 = -qJ(3) - pkin(6);
t31 = t36 * pkin(2) + pkin(1);
t1 = [(-g(1) * (-t35 * t31 - t37 * t33) - g(2) * (t37 * t31 - t35 * t33)) * MDP(12) + (MDP(3) - MDP(11)) * t39 + (-t34 * MDP(10) + MDP(18) * t29 - MDP(19) * t28 + MDP(25) * t27 - MDP(26) * t26 + t36 * MDP(9) + MDP(2)) * t24; (g(3) * t34 + t39 * t36) * MDP(10) + t40 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t36 + t39 * t34); -t24 * MDP(12); t40; t41;];
taug = t1;
