% Calculate Gravitation load on the joints for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:23
% EndTime: 2019-12-05 17:03:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (95->31), mult. (139->47), div. (0->0), fcn. (126->8), ass. (0->21)
t36 = qJ(3) + qJ(4);
t34 = sin(t36);
t52 = g(2) * t34;
t37 = sin(qJ(5));
t39 = sin(qJ(2));
t50 = t39 * t37;
t40 = cos(qJ(5));
t49 = t39 * t40;
t42 = cos(qJ(2));
t48 = t42 * t37;
t47 = t42 * t40;
t35 = cos(t36);
t45 = g(1) * t42 + g(3) * t39;
t46 = (t45 * t35 - t52) * MDP(18) + (t40 * MDP(24) - t37 * MDP(25) + MDP(17)) * (g(2) * t35 + t45 * t34);
t41 = cos(qJ(3));
t38 = sin(qJ(3));
t33 = t35 * t47 + t50;
t32 = -t35 * t48 + t49;
t31 = -t35 * t49 + t48;
t30 = t35 * t50 + t47;
t1 = [-g(3) * MDP(1); t45 * MDP(4) + (-g(1) * t31 - g(3) * t33) * MDP(24) + (-g(1) * t30 - g(3) * t32) * MDP(25) + (t41 * MDP(10) - t38 * MDP(11) + MDP(17) * t35 - MDP(18) * t34 + MDP(3)) * (g(1) * t39 - g(3) * t42); (g(2) * t41 + t45 * t38) * MDP(10) + (-g(2) * t38 + t45 * t41) * MDP(11) + t46; t46; (-g(1) * t32 + g(3) * t30 - t37 * t52) * MDP(24) + (g(1) * t33 - g(3) * t31 - t40 * t52) * MDP(25);];
taug = t1;
