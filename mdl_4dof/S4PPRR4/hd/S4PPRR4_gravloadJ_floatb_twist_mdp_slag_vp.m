% Calculate Gravitation load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (40->18), mult. (57->31), div. (0->0), fcn. (52->6), ass. (0->14)
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t31 = t22 * MDP(11) - t21 * MDP(12) + MDP(4);
t18 = pkin(7) + qJ(3);
t16 = sin(t18);
t30 = g(3) * t16;
t19 = sin(pkin(6));
t29 = t19 * t21;
t28 = t19 * t22;
t20 = cos(pkin(6));
t27 = t20 * t21;
t26 = t20 * t22;
t17 = cos(t18);
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t19 + g(2) * t20) * MDP(2); (MDP(5) * t16 - t31 * t17) * g(3) + (MDP(5) * t17 + t31 * t16) * (g(1) * t20 + g(2) * t19); (-g(1) * (-t17 * t27 + t28) - g(2) * (-t17 * t29 - t26) + t21 * t30) * MDP(11) + (-g(1) * (-t17 * t26 - t29) - g(2) * (-t17 * t28 + t27) + t22 * t30) * MDP(12);];
taug = t1;
