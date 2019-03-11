% Calculate Gravitation load on the joints for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:56
% EndTime: 2019-03-08 18:17:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (55->20), mult. (46->21), div. (0->0), fcn. (28->2), ass. (0->9)
t19 = pkin(5) + qJ(2);
t17 = sin(t19);
t18 = cos(t19);
t21 = t18 * pkin(2) + t17 * qJ(3);
t20 = -MDP(10) - MDP(7);
t14 = t18 * qJ(3);
t12 = g(1) * t18 + g(2) * t17;
t11 = g(1) * t17 - g(2) * t18;
t1 = [(-MDP(1) + t20) * g(3); (-g(1) * (-t17 * pkin(2) + t14) - g(2) * t21) * MDP(7) + (-g(1) * (t14 + (-pkin(2) - qJ(4)) * t17) - g(2) * (t18 * qJ(4) + t21)) * MDP(10) + (MDP(4) - MDP(6) - MDP(8)) * t12 + (MDP(3) - MDP(5) + MDP(9)) * t11; t20 * t11; -t12 * MDP(10);];
taug  = t1;
