% Calculate Gravitation load on the joints for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:50
% EndTime: 2019-03-08 18:31:50
% DurationCPUTime: 0.03s
% Computational Cost: add. (76->15), mult. (42->20), div. (0->0), fcn. (26->6), ass. (0->11)
t23 = qJ(1) + pkin(7) + qJ(3);
t22 = qJ(4) + t23;
t18 = sin(t22);
t19 = cos(t22);
t28 = (g(1) * t18 - g(2) * t19) * MDP(9) + (g(1) * t19 + g(2) * t18) * MDP(10);
t20 = sin(t23);
t21 = cos(t23);
t27 = (g(1) * t20 - g(2) * t21) * MDP(6) + (g(1) * t21 + g(2) * t20) * MDP(7) + t28;
t25 = cos(qJ(1));
t24 = sin(qJ(1));
t1 = [(g(1) * t25 + g(2) * t24) * MDP(3) + t27 + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t24 - g(2) * t25); -g(3) * MDP(4); t27; t28;];
taug  = t1;
