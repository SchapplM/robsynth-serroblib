% Calculate Gravitation load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:18
% EndTime: 2019-07-18 13:30:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (103->19), mult. (79->27), div. (0->0), fcn. (58->8), ass. (0->14)
t28 = qJ(2) + qJ(3);
t27 = qJ(4) + t28;
t23 = sin(t27);
t24 = cos(t27);
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t35 = g(1) * t24 + g(2) * t23;
t36 = t35 * MDP(10) + (t31 * MDP(16) - t29 * MDP(17) + MDP(9)) * (g(1) * t23 - g(2) * t24);
t25 = sin(t28);
t26 = cos(t28);
t33 = (g(1) * t25 - g(2) * t26) * MDP(6) + (g(1) * t26 + g(2) * t25) * MDP(7) + t36;
t32 = cos(qJ(2));
t30 = sin(qJ(2));
t1 = [-g(3) * MDP(1); (g(1) * t30 - g(2) * t32) * MDP(3) + (g(1) * t32 + g(2) * t30) * MDP(4) + t33; t33; t36; (-g(3) * t31 + t35 * t29) * MDP(16) + (g(3) * t29 + t35 * t31) * MDP(17);];
taug  = t1;
