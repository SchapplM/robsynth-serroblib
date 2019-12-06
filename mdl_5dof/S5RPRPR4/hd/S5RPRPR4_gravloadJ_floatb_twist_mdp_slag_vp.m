% Calculate Gravitation load on the joints for
% S5RPRPR4
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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:58
% EndTime: 2019-12-05 17:53:58
% DurationCPUTime: 0.09s
% Computational Cost: add. (108->31), mult. (95->42), div. (0->0), fcn. (69->8), ass. (0->16)
t26 = qJ(3) + pkin(9) + qJ(5);
t21 = sin(t26);
t22 = cos(t26);
t27 = qJ(1) + pkin(8);
t24 = sin(t27);
t25 = cos(t27);
t35 = g(2) * t24 - g(3) * t25;
t37 = (-g(1) * t22 - t35 * t21) * MDP(19) + (g(1) * t21 - t35 * t22) * MDP(20);
t36 = g(2) * t25 + g(3) * t24;
t32 = cos(qJ(1));
t31 = cos(qJ(3));
t30 = sin(qJ(1));
t29 = sin(qJ(3));
t28 = -qJ(4) - pkin(6);
t23 = t31 * pkin(3) + pkin(2);
t1 = [(-g(2) * t30 + g(3) * t32) * MDP(3) + t35 * MDP(12) + (-g(2) * (-t32 * pkin(1) - t25 * t23 + t24 * t28) - g(3) * (-t30 * pkin(1) - t24 * t23 - t25 * t28)) * MDP(13) + (pkin(1) * MDP(4) + MDP(2)) * (g(2) * t32 + g(3) * t30) + (t31 * MDP(10) - t29 * MDP(11) + MDP(19) * t22 - MDP(20) * t21) * t36; (-MDP(13) - MDP(4)) * g(1); (g(1) * t29 - t35 * t31) * MDP(11) + t37 + (MDP(13) * pkin(3) + MDP(10)) * (-g(1) * t31 - t35 * t29); -t36 * MDP(13); t37;];
taug = t1;
