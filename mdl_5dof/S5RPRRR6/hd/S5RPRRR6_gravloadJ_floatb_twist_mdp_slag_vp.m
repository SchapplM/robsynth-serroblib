% Calculate Gravitation load on the joints for
% S5RPRRR6
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (141->35), mult. (144->55), div. (0->0), fcn. (128->10), ass. (0->22)
t41 = qJ(3) + qJ(4);
t38 = sin(t41);
t56 = g(3) * t38;
t39 = cos(t41);
t42 = sin(qJ(5));
t54 = t39 * t42;
t45 = cos(qJ(5));
t53 = t39 * t45;
t40 = qJ(1) + pkin(9);
t36 = sin(t40);
t37 = cos(t40);
t51 = g(1) * t37 + g(2) * t36;
t52 = (t51 * t39 + t56) * MDP(18) + (t45 * MDP(24) - t42 * MDP(25) + MDP(17)) * (-g(3) * t39 + t51 * t38);
t47 = cos(qJ(1));
t46 = cos(qJ(3));
t44 = sin(qJ(1));
t43 = sin(qJ(3));
t35 = t36 * t42 + t37 * t53;
t34 = t36 * t45 - t37 * t54;
t33 = -t36 * t53 + t37 * t42;
t32 = t36 * t54 + t37 * t45;
t1 = [(g(1) * t47 + g(2) * t44) * MDP(3) + (-g(1) * t33 - g(2) * t35) * MDP(24) + (-g(1) * t32 - g(2) * t34) * MDP(25) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t44 - g(2) * t47) + (t46 * MDP(10) - t43 * MDP(11) + MDP(17) * t39 - MDP(18) * t38) * (g(1) * t36 - g(2) * t37); -g(3) * MDP(4); (-g(3) * t46 + t51 * t43) * MDP(10) + (g(3) * t43 + t51 * t46) * MDP(11) + t52; t52; (-g(1) * t34 + g(2) * t32 + t42 * t56) * MDP(24) + (g(1) * t35 - g(2) * t33 + t45 * t56) * MDP(25);];
taug = t1;
