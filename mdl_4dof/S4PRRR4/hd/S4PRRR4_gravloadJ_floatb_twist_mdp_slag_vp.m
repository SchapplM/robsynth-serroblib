% Calculate Gravitation load on the joints for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:39
% DurationCPUTime: 0.05s
% Computational Cost: add. (68->17), mult. (63->23), div. (0->0), fcn. (50->6), ass. (0->11)
t20 = qJ(3) + qJ(4);
t17 = sin(t20);
t18 = cos(t20);
t19 = pkin(7) + qJ(2);
t15 = sin(t19);
t16 = cos(t19);
t24 = g(1) * t16 + g(2) * t15;
t25 = (-g(3) * t18 + t24 * t17) * MDP(17) + (g(3) * t17 + t24 * t18) * MDP(18);
t22 = cos(qJ(3));
t21 = sin(qJ(3));
t1 = [-g(3) * MDP(1); t24 * MDP(4) + (t22 * MDP(10) - t21 * MDP(11) + MDP(17) * t18 - MDP(18) * t17 + MDP(3)) * (g(1) * t15 - g(2) * t16); (-g(3) * t22 + t24 * t21) * MDP(10) + (g(3) * t21 + t24 * t22) * MDP(11) + t25; t25;];
taug = t1;
