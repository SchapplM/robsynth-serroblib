% Calculate Gravitation load on the joints for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:08
% EndTime: 2019-03-08 18:21:08
% DurationCPUTime: 0.04s
% Computational Cost: add. (58->18), mult. (52->25), div. (0->0), fcn. (46->4), ass. (0->10)
t24 = pkin(6) + qJ(2);
t22 = sin(t24);
t23 = cos(t24);
t27 = sin(qJ(4));
t28 = cos(qJ(4));
t15 = -t22 * t27 - t23 * t28;
t16 = -t22 * t28 + t23 * t27;
t29 = -(g(1) * t15 + g(2) * t16) * MDP(10) + (g(1) * t16 - g(2) * t15) * MDP(9);
t17 = g(1) * t22 - g(2) * t23;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(1) * (-t22 * pkin(2) + t23 * qJ(3)) - g(2) * (t23 * pkin(2) + t22 * qJ(3))) * MDP(7) + (MDP(4) - MDP(6)) * (g(1) * t23 + g(2) * t22) + (MDP(3) + MDP(5)) * t17 - t29; -t17 * MDP(7); t29;];
taug  = t1;
