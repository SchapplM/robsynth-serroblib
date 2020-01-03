% Calculate Gravitation load on the joints for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:49
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->26), mult. (87->38), div. (0->0), fcn. (80->6), ass. (0->13)
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t36 = t32 * pkin(1) + t30 * qJ(2);
t35 = cos(pkin(6));
t28 = sin(pkin(6));
t19 = t32 * t28 - t30 * t35;
t20 = t30 * t28 + t32 * t35;
t34 = g(1) * t20 - g(2) * t19;
t31 = cos(qJ(4));
t29 = sin(qJ(4));
t25 = t32 * qJ(2);
t21 = g(1) * t30 - g(2) * t32;
t1 = [(-g(1) * (-t30 * pkin(1) + t25) - g(2) * t36) * MDP(6) - t34 * MDP(8) + (-g(1) * (t25 + (-pkin(1) - pkin(2)) * t30) - g(2) * (t32 * pkin(2) + t36)) * MDP(9) + (MDP(3) - MDP(5)) * (g(1) * t32 + g(2) * t30) + (MDP(2) + MDP(4)) * t21 + (-t31 * MDP(15) + t29 * MDP(16) - MDP(7)) * (g(1) * t19 + g(2) * t20); (-MDP(6) - MDP(9)) * t21; g(3) * MDP(9); (g(3) * t31 + t34 * t29) * MDP(15) + (-g(3) * t29 + t34 * t31) * MDP(16);];
taug = t1;
