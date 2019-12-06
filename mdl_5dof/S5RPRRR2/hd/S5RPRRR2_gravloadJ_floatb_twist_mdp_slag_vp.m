% Calculate Gravitation load on the joints for
% S5RPRRR2
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:14
% EndTime: 2019-12-05 18:12:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (167->31), mult. (131->40), div. (0->0), fcn. (104->10), ass. (0->16)
t34 = pkin(9) + qJ(3);
t33 = qJ(4) + t34;
t30 = qJ(5) + t33;
t26 = sin(t30);
t27 = cos(t30);
t37 = sin(qJ(1));
t38 = cos(qJ(1));
t39 = g(1) * t38 + g(2) * t37;
t41 = (-g(3) * t27 + t39 * t26) * MDP(27) + (g(3) * t26 + t39 * t27) * MDP(28);
t28 = sin(t33);
t29 = cos(t33);
t40 = (-g(3) * t29 + t39 * t28) * MDP(20) + (g(3) * t28 + t39 * t29) * MDP(21) + t41;
t24 = g(1) * t37 - g(2) * t38;
t32 = cos(t34);
t31 = sin(t34);
t1 = [(-g(1) * (-t37 * pkin(1) + t38 * qJ(2)) - g(2) * (t38 * pkin(1) + t37 * qJ(2))) * MDP(7) + (MDP(3) - MDP(6)) * t39 + (t32 * MDP(13) - t31 * MDP(14) + MDP(20) * t29 - MDP(21) * t28 + MDP(27) * t27 - MDP(28) * t26 + MDP(4) * cos(pkin(9)) - MDP(5) * sin(pkin(9)) + MDP(2)) * t24; -t24 * MDP(7); (-g(3) * t32 + t39 * t31) * MDP(13) + (g(3) * t31 + t39 * t32) * MDP(14) + t40; t40; t41;];
taug = t1;
