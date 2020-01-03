% Calculate Gravitation load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:42
% EndTime: 2020-01-03 11:47:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (112->34), mult. (109->49), div. (0->0), fcn. (79->8), ass. (0->18)
t34 = qJ(3) + qJ(4);
t29 = sin(t34);
t30 = cos(t34);
t33 = qJ(1) + pkin(8);
t27 = sin(t33);
t28 = cos(t33);
t41 = g(2) * t27 - g(3) * t28;
t39 = -g(1) * t30 + t41 * t29;
t44 = t39 * MDP(17) + (g(1) * t29 + t41 * t30) * MDP(18);
t37 = cos(qJ(3));
t43 = t37 * pkin(3) + pkin(4) * t30;
t42 = g(2) * t28 + g(3) * t27;
t38 = cos(qJ(1));
t36 = sin(qJ(1));
t35 = sin(qJ(3));
t32 = -qJ(5) - pkin(7) - pkin(6);
t24 = pkin(2) + t43;
t1 = [(g(2) * t36 - g(3) * t38) * MDP(3) - t41 * MDP(19) + (-g(2) * (t38 * pkin(1) + t28 * t24 - t27 * t32) - g(3) * (t36 * pkin(1) + t27 * t24 + t28 * t32)) * MDP(20) + (pkin(1) * MDP(4) + MDP(2)) * (-g(2) * t38 - g(3) * t36) + (-t37 * MDP(10) + t35 * MDP(11) - MDP(17) * t30 + MDP(18) * t29) * t42; (-MDP(20) - MDP(4)) * g(1); (-g(1) * t37 + t41 * t35) * MDP(10) + (g(1) * t35 + t41 * t37) * MDP(11) + (-g(1) * t43 - t41 * (-t35 * pkin(3) - pkin(4) * t29)) * MDP(20) + t44; t39 * MDP(20) * pkin(4) + t44; t42 * MDP(20);];
taug = t1;
