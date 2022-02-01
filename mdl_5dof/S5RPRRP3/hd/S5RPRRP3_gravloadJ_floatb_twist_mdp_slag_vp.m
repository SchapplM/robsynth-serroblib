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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:20
% EndTime: 2022-01-23 09:30:20
% DurationCPUTime: 0.14s
% Computational Cost: add. (156->37), mult. (143->49), div. (0->0), fcn. (107->8), ass. (0->20)
t53 = MDP(17) + MDP(19);
t52 = MDP(18) + MDP(20);
t42 = qJ(3) + qJ(4);
t38 = cos(t42);
t45 = cos(qJ(3));
t51 = t45 * pkin(3) + pkin(4) * t38;
t37 = sin(t42);
t41 = qJ(1) + pkin(8);
t35 = sin(t41);
t36 = cos(t41);
t49 = g(1) * t36 + g(2) * t35;
t28 = -g(3) * t38 + t49 * t37;
t50 = t52 * (g(3) * t37 + t49 * t38) + t53 * t28;
t48 = g(1) * t35 - g(2) * t36;
t46 = cos(qJ(1));
t44 = sin(qJ(1));
t43 = sin(qJ(3));
t40 = -qJ(5) - pkin(7) - pkin(6);
t32 = pkin(2) + t51;
t1 = [(g(1) * t46 + g(2) * t44) * MDP(3) - t49 * MDP(21) + (-g(1) * (-t44 * pkin(1) - t35 * t32 - t36 * t40) - g(2) * (t46 * pkin(1) + t36 * t32 - t35 * t40)) * MDP(22) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t44 - g(2) * t46) + (t45 * MDP(10) - t43 * MDP(11) - t52 * t37 + t53 * t38) * t48; (-MDP(22) - MDP(4)) * g(3); (-g(3) * t45 + t49 * t43) * MDP(10) + (g(3) * t43 + t49 * t45) * MDP(11) + (-g(3) * t51 - t49 * (-t43 * pkin(3) - pkin(4) * t37)) * MDP(22) + t50; t28 * MDP(22) * pkin(4) + t50; -t48 * MDP(22);];
taug = t1;
