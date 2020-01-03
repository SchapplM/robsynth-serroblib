% Calculate Gravitation load on the joints for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (94->30), mult. (138->46), div. (0->0), fcn. (126->8), ass. (0->21)
t36 = qJ(2) + qJ(3);
t34 = sin(t36);
t52 = g(3) * t34;
t37 = sin(qJ(4));
t39 = sin(qJ(1));
t50 = t39 * t37;
t40 = cos(qJ(4));
t49 = t39 * t40;
t42 = cos(qJ(1));
t48 = t42 * t37;
t47 = t42 * t40;
t35 = cos(t36);
t45 = g(1) * t42 + g(2) * t39;
t46 = (t45 * t35 + t52) * MDP(17) + (t40 * MDP(23) - t37 * MDP(24) + MDP(16)) * (-g(3) * t35 + t45 * t34);
t41 = cos(qJ(2));
t38 = sin(qJ(2));
t33 = t35 * t47 + t50;
t32 = -t35 * t48 + t49;
t31 = -t35 * t49 + t48;
t30 = t35 * t50 + t47;
t1 = [t45 * MDP(3) + (-g(1) * t31 - g(2) * t33) * MDP(23) + (-g(1) * t30 - g(2) * t32) * MDP(24) + (-t38 * MDP(10) + MDP(16) * t35 - MDP(17) * t34 + t41 * MDP(9) + MDP(2)) * (g(1) * t39 - g(2) * t42); (-g(3) * t41 + t45 * t38) * MDP(9) + (g(3) * t38 + t45 * t41) * MDP(10) + t46; t46; (-g(1) * t32 + g(2) * t30 + t37 * t52) * MDP(23) + (g(1) * t33 - g(2) * t31 + t40 * t52) * MDP(24);];
taug = t1;
