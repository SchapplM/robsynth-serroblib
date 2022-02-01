% Calculate Gravitation load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:21
% DurationCPUTime: 0.09s
% Computational Cost: add. (141->24), mult. (91->36), div. (0->0), fcn. (64->8), ass. (0->15)
t28 = qJ(1) + qJ(2);
t25 = pkin(9) + qJ(4) + t28;
t23 = sin(t25);
t24 = cos(t25);
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t36 = g(1) * t24 + g(2) * t23;
t37 = t36 * MDP(10) + (t31 * MDP(16) - t29 * MDP(17) + MDP(9)) * (g(1) * t23 - g(2) * t24);
t26 = sin(t28);
t27 = cos(t28);
t34 = g(1) * t26 - g(2) * t27;
t33 = t34 * MDP(5) + (g(1) * t27 + g(2) * t26) * MDP(6) + t37;
t32 = cos(qJ(1));
t30 = sin(qJ(1));
t1 = [(g(1) * t30 - g(2) * t32) * MDP(2) + (g(1) * t32 + g(2) * t30) * MDP(3) + (-g(1) * (-t30 * pkin(1) - pkin(2) * t26) - g(2) * (t32 * pkin(1) + pkin(2) * t27)) * MDP(7) + t33; t34 * MDP(7) * pkin(2) + t33; -g(3) * MDP(7); t37; (-g(3) * t31 + t36 * t29) * MDP(16) + (g(3) * t29 + t36 * t31) * MDP(17);];
taug = t1;
