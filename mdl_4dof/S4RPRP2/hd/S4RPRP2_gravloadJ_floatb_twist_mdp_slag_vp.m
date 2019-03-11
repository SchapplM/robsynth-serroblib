% Calculate Gravitation load on the joints for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:53
% EndTime: 2019-03-08 18:30:53
% DurationCPUTime: 0.07s
% Computational Cost: add. (43->26), mult. (80->35), div. (0->0), fcn. (68->4), ass. (0->15)
t36 = cos(qJ(3));
t37 = cos(qJ(1));
t34 = sin(qJ(3));
t35 = sin(qJ(1));
t42 = t35 * t34;
t22 = -t37 * t36 - t42;
t41 = t37 * t34;
t23 = -t35 * t36 + t41;
t43 = (g(1) * t22 + g(2) * t23) * MDP(9);
t40 = t37 * pkin(1) + t35 * qJ(2);
t39 = g(1) * t23 - g(2) * t22;
t31 = t37 * qJ(2);
t29 = t36 * pkin(3) + pkin(2);
t24 = g(1) * t35 - g(2) * t37;
t1 = [(-g(1) * (-t35 * pkin(1) + t31) - g(2) * t40) * MDP(6) - t39 * MDP(8) + t43 + (-g(1) * (pkin(3) * t41 + t31 + (-pkin(1) - t29) * t35) - g(2) * (pkin(3) * t42 + t37 * t29 + t40)) * MDP(10) + (MDP(3) - MDP(5)) * (g(1) * t37 + g(2) * t35) + (MDP(2) + MDP(4)) * t24; (-MDP(10) - MDP(6)) * t24; -t43 + (MDP(10) * pkin(3) + MDP(8)) * t39; g(3) * MDP(10);];
taug  = t1;
