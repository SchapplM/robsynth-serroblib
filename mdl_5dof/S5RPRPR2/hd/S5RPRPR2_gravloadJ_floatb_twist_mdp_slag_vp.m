% Calculate Gravitation load on the joints for
% S5RPRPR2
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:07
% DurationCPUTime: 0.07s
% Computational Cost: add. (139->30), mult. (94->38), div. (0->0), fcn. (66->9), ass. (0->15)
t39 = qJ(1) + pkin(8);
t37 = qJ(3) + t39;
t33 = sin(t37);
t34 = cos(t37);
t46 = t34 * pkin(3) + t33 * qJ(4);
t45 = -t33 * pkin(3) + t34 * qJ(4);
t29 = g(1) * t34 + g(2) * t33;
t28 = g(1) * t33 - g(2) * t34;
t38 = pkin(9) + qJ(5);
t35 = sin(t38);
t36 = cos(t38);
t43 = (MDP(7) - MDP(9)) * t29 + (t36 * MDP(16) - t35 * MDP(17) + MDP(8) * cos(pkin(9)) + MDP(6)) * t28;
t42 = cos(qJ(1));
t41 = sin(qJ(1));
t1 = [(g(1) * t42 + g(2) * t41) * MDP(3) + (-g(1) * (-pkin(2) * sin(t39) - t41 * pkin(1) + t45) - g(2) * (pkin(2) * cos(t39) + t42 * pkin(1) + t46)) * MDP(10) + t43 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t41 - g(2) * t42); (-MDP(10) - MDP(4)) * g(3); (-g(1) * t45 - g(2) * t46) * MDP(10) + t43; -t28 * MDP(10); (-g(3) * t36 + t29 * t35) * MDP(16) + (g(3) * t35 + t29 * t36) * MDP(17);];
taug = t1;
