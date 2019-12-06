% Calculate Gravitation load on the joints for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:41
% EndTime: 2019-12-05 16:19:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (108->26), mult. (88->33), div. (0->0), fcn. (65->6), ass. (0->14)
t25 = qJ(3) + pkin(9) + qJ(5);
t20 = sin(t25);
t21 = cos(t25);
t26 = pkin(8) + qJ(2);
t23 = sin(t26);
t24 = cos(t26);
t31 = g(1) * t24 + g(2) * t23;
t32 = (-g(3) * t21 + t31 * t20) * MDP(19) + (g(3) * t20 + t31 * t21) * MDP(20);
t18 = g(1) * t23 - g(2) * t24;
t29 = cos(qJ(3));
t28 = sin(qJ(3));
t27 = -qJ(4) - pkin(6);
t22 = t29 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(1) * (-t23 * t22 - t24 * t27) - g(2) * (t24 * t22 - t23 * t27)) * MDP(13) + (MDP(4) - MDP(12)) * t31 + (t29 * MDP(10) - t28 * MDP(11) + MDP(19) * t21 - MDP(20) * t20 + MDP(3)) * t18; (g(3) * t28 + t31 * t29) * MDP(11) + t32 + (MDP(13) * pkin(3) + MDP(10)) * (-g(3) * t29 + t31 * t28); -t18 * MDP(13); t32;];
taug = t1;
