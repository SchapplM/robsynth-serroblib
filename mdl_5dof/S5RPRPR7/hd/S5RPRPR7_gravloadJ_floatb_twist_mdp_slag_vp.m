% Calculate Gravitation load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:51
% DurationCPUTime: 0.15s
% Computational Cost: add. (108->41), mult. (119->62), div. (0->0), fcn. (101->10), ass. (0->26)
t39 = qJ(3) + pkin(9);
t35 = sin(t39);
t57 = g(3) * t35;
t40 = qJ(1) + pkin(8);
t36 = sin(t40);
t42 = sin(qJ(5));
t55 = t36 * t42;
t45 = cos(qJ(5));
t54 = t36 * t45;
t38 = cos(t40);
t53 = t38 * t42;
t52 = t38 * t45;
t51 = g(1) * t38 + g(2) * t36;
t50 = g(1) * t36 - g(2) * t38;
t47 = cos(qJ(1));
t46 = cos(qJ(3));
t44 = sin(qJ(1));
t43 = sin(qJ(3));
t41 = -qJ(4) - pkin(6);
t37 = cos(t39);
t34 = t46 * pkin(3) + pkin(2);
t33 = t37 * t52 + t55;
t32 = -t37 * t53 + t54;
t31 = -t37 * t54 + t53;
t30 = t37 * t55 + t52;
t1 = [(g(1) * t47 + g(2) * t44) * MDP(3) - t51 * MDP(12) + (-g(1) * (-t44 * pkin(1) - t36 * t34 - t38 * t41) - g(2) * (t47 * pkin(1) + t38 * t34 - t36 * t41)) * MDP(13) + (-g(1) * t31 - g(2) * t33) * MDP(19) + (-g(1) * t30 - g(2) * t32) * MDP(20) + (t46 * MDP(10) - t43 * MDP(11)) * t50 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t44 - g(2) * t47); (-MDP(13) - MDP(4)) * g(3); (g(3) * t43 + t51 * t46) * MDP(11) + (MDP(13) * pkin(3) + MDP(10)) * (-g(3) * t46 + t51 * t43) + (MDP(19) * t45 - MDP(20) * t42) * (-g(3) * t37 + t51 * t35); -t50 * MDP(13); (-g(1) * t32 + g(2) * t30 + t42 * t57) * MDP(19) + (g(1) * t33 - g(2) * t31 + t45 * t57) * MDP(20);];
taug = t1;
