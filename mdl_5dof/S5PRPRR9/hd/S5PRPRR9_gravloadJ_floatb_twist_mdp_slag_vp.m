% Calculate Gravitation load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (107->23), mult. (108->33), div. (0->0), fcn. (108->6), ass. (0->13)
t32 = pkin(8) + qJ(2);
t30 = sin(t32);
t31 = cos(t32);
t35 = sin(qJ(4));
t36 = cos(qJ(4));
t19 = -t30 * t35 - t31 * t36;
t20 = -t30 * t36 + t31 * t35;
t26 = sin(qJ(5));
t27 = cos(qJ(5));
t28 = g(1) * t19 + g(2) * t20;
t39 = (t27 * MDP(16) - t26 * MDP(17) + MDP(9)) * (g(1) * t20 - g(2) * t19) - t28 * MDP(10);
t21 = g(1) * t30 - g(2) * t31;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(1) * (-t30 * pkin(2) + t31 * qJ(3)) - g(2) * (t31 * pkin(2) + t30 * qJ(3))) * MDP(7) + (MDP(4) - MDP(6)) * (g(1) * t31 + g(2) * t30) + (MDP(3) + MDP(5)) * t21 - t39; -t21 * MDP(7); t39; (g(3) * t27 - t28 * t26) * MDP(16) + (-g(3) * t26 - t28 * t27) * MDP(17);];
taug = t1;
