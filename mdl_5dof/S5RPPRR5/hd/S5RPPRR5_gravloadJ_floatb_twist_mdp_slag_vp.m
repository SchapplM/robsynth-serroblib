% Calculate Gravitation load on the joints for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:39
% EndTime: 2019-12-31 17:56:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (107->28), mult. (115->42), div. (0->0), fcn. (112->8), ass. (0->15)
t36 = qJ(1) + pkin(8);
t34 = sin(t36);
t35 = cos(t36);
t39 = sin(qJ(4));
t40 = cos(qJ(4));
t21 = -t34 * t39 - t35 * t40;
t22 = -t34 * t40 + t35 * t39;
t27 = sin(qJ(5));
t29 = cos(qJ(5));
t32 = g(1) * t21 + g(2) * t22;
t43 = (t29 * MDP(16) - t27 * MDP(17) + MDP(9)) * (g(1) * t22 - g(2) * t21) - t32 * MDP(10);
t30 = cos(qJ(1));
t28 = sin(qJ(1));
t23 = g(1) * t34 - g(2) * t35;
t1 = [(g(1) * t30 + g(2) * t28) * MDP(3) + t23 * MDP(5) + (-g(1) * t35 - g(2) * t34) * MDP(6) + (-g(1) * (-t28 * pkin(1) - t34 * pkin(2) + t35 * qJ(3)) - g(2) * (t30 * pkin(1) + t35 * pkin(2) + t34 * qJ(3))) * MDP(7) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t28 - g(2) * t30) - t43; (-MDP(4) - MDP(7)) * g(3); -t23 * MDP(7); t43; (g(3) * t29 - t32 * t27) * MDP(16) + (-g(3) * t27 - t32 * t29) * MDP(17);];
taug = t1;
