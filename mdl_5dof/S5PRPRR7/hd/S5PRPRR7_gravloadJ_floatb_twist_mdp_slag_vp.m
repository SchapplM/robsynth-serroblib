% Calculate Gravitation load on the joints for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:46
% EndTime: 2019-12-05 16:00:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (86->36), mult. (148->57), div. (0->0), fcn. (139->8), ass. (0->18)
t34 = sin(pkin(8));
t35 = cos(pkin(8));
t42 = g(1) * t35 + g(2) * t34;
t39 = cos(qJ(2));
t48 = g(3) * t39;
t37 = sin(qJ(2));
t47 = t34 * t37;
t46 = t35 * t37;
t36 = sin(qJ(4));
t45 = t36 * t37;
t38 = cos(qJ(4));
t44 = t37 * t38;
t33 = qJ(4) + qJ(5);
t30 = sin(t33);
t31 = cos(t33);
t43 = (-g(1) * (-t30 * t34 + t31 * t46) - g(2) * (t30 * t35 + t31 * t47) + t31 * t48) * MDP(20) + (-g(1) * (-t30 * t46 - t31 * t34) - g(2) * (-t30 * t47 + t31 * t35) - t30 * t48) * MDP(21);
t28 = t37 * t42 - t48;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(3) * (pkin(2) * t39 + qJ(3) * t37) + t42 * (pkin(2) * t37 - qJ(3) * t39)) * MDP(7) + (MDP(3) - MDP(5)) * t28 + (-MDP(13) * t36 - MDP(14) * t38 - MDP(20) * t30 - MDP(21) * t31 + MDP(4) - MDP(6)) * (g(3) * t37 + t39 * t42); -t28 * MDP(7); (-g(1) * (-t34 * t36 + t35 * t44) - g(2) * (t34 * t44 + t35 * t36) + t38 * t48) * MDP(13) + (-g(1) * (-t34 * t38 - t35 * t45) - g(2) * (-t34 * t45 + t35 * t38) - t36 * t48) * MDP(14) + t43; t43;];
taug = t1;
