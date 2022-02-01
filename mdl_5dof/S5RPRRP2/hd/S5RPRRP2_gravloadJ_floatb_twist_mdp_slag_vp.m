% Calculate Gravitation load on the joints for
% S5RPRRP2
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:13
% EndTime: 2022-01-23 09:28:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (166->33), mult. (129->39), div. (0->0), fcn. (93->8), ass. (0->18)
t51 = MDP(13) + MDP(15);
t50 = MDP(14) + MDP(16);
t39 = qJ(1) + pkin(8);
t38 = qJ(3) + t39;
t35 = sin(t38);
t36 = cos(t38);
t43 = cos(qJ(4));
t37 = t43 * pkin(4) + pkin(3);
t40 = -qJ(5) - pkin(7);
t49 = -t35 * t40 + t36 * t37;
t33 = g(1) * t36 + g(2) * t35;
t32 = g(1) * t35 - g(2) * t36;
t47 = -t35 * t37 - t36 * t40;
t41 = sin(qJ(4));
t46 = (-MDP(17) + MDP(7)) * t33 + (-t50 * t41 + t51 * t43 + MDP(6)) * t32;
t44 = cos(qJ(1));
t42 = sin(qJ(1));
t1 = [(g(1) * t44 + g(2) * t42) * MDP(3) + (-g(1) * (-pkin(2) * sin(t39) - t42 * pkin(1) + t47) - g(2) * (pkin(2) * cos(t39) + t44 * pkin(1) + t49)) * MDP(18) + t46 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t42 - g(2) * t44); (-MDP(18) - MDP(4)) * g(3); (-g(1) * t47 - g(2) * t49) * MDP(18) + t46; t50 * (g(3) * t41 + t33 * t43) + (MDP(18) * pkin(4) + t51) * (-g(3) * t43 + t33 * t41); -t32 * MDP(18);];
taug = t1;
