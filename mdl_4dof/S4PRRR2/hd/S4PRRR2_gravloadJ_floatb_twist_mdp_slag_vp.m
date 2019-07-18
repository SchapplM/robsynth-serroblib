% Calculate Gravitation load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:24
% EndTime: 2019-07-18 13:27:24
% DurationCPUTime: 0.03s
% Computational Cost: add. (53->13), mult. (37->19), div. (0->0), fcn. (24->6), ass. (0->11)
t22 = qJ(2) + qJ(3);
t21 = qJ(4) + t22;
t17 = sin(t21);
t18 = cos(t21);
t26 = (g(1) * t18 + g(3) * t17) * MDP(9) + (-g(1) * t17 + g(3) * t18) * MDP(10);
t19 = sin(t22);
t20 = cos(t22);
t25 = (g(1) * t20 + g(3) * t19) * MDP(6) + (-g(1) * t19 + g(3) * t20) * MDP(7) + t26;
t24 = cos(qJ(2));
t23 = sin(qJ(2));
t1 = [g(2) * MDP(1); (g(1) * t24 + g(3) * t23) * MDP(3) + (-g(1) * t23 + g(3) * t24) * MDP(4) + t25; t25; t26;];
taug  = t1;
