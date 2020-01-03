% Calculate Gravitation load on the joints for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (65->34), mult. (110->52), div. (0->0), fcn. (97->8), ass. (0->23)
t35 = qJ(2) + pkin(7);
t33 = sin(t35);
t50 = g(3) * t33;
t37 = sin(qJ(4));
t39 = sin(qJ(1));
t48 = t39 * t37;
t40 = cos(qJ(4));
t47 = t39 * t40;
t42 = cos(qJ(1));
t46 = t42 * t37;
t45 = t42 * t40;
t44 = g(1) * t42 + g(2) * t39;
t30 = g(1) * t39 - g(2) * t42;
t41 = cos(qJ(2));
t38 = sin(qJ(2));
t36 = -qJ(3) - pkin(5);
t34 = cos(t35);
t32 = t41 * pkin(2) + pkin(1);
t29 = t34 * t45 + t48;
t28 = -t34 * t46 + t47;
t27 = -t34 * t47 + t46;
t26 = t34 * t48 + t45;
t1 = [(-g(1) * (-t39 * t32 - t42 * t36) - g(2) * (t42 * t32 - t39 * t36)) * MDP(12) + (-g(1) * t27 - g(2) * t29) * MDP(18) + (-g(1) * t26 - g(2) * t28) * MDP(19) + (MDP(3) - MDP(11)) * t44 + (-t38 * MDP(10) + t41 * MDP(9) + MDP(2)) * t30; (g(3) * t38 + t44 * t41) * MDP(10) + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t41 + t44 * t38) + (MDP(18) * t40 - MDP(19) * t37) * (-g(3) * t34 + t44 * t33); -t30 * MDP(12); (-g(1) * t28 + g(2) * t26 + t37 * t50) * MDP(18) + (g(1) * t29 - g(2) * t27 + t40 * t50) * MDP(19);];
taug = t1;
