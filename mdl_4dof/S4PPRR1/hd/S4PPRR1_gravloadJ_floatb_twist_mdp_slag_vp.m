% Calculate Gravitation load on the joints for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:49
% EndTime: 2019-03-08 18:15:49
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->14), mult. (47->24), div. (0->0), fcn. (50->6), ass. (0->13)
t22 = qJ(3) + qJ(4);
t20 = sin(t22);
t21 = cos(t22);
t23 = sin(pkin(6));
t24 = cos(pkin(6));
t16 = -t23 * t20 - t24 * t21;
t17 = t24 * t20 - t23 * t21;
t27 = (g(1) * t17 - g(2) * t16) * MDP(7) + (-g(1) * t16 - g(2) * t17) * MDP(8);
t26 = cos(qJ(3));
t25 = sin(qJ(3));
t19 = -t23 * t26 + t24 * t25;
t18 = -t23 * t25 - t24 * t26;
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t23 + g(2) * t24) * MDP(2); (g(1) * t19 - g(2) * t18) * MDP(4) + (-g(1) * t18 - g(2) * t19) * MDP(5) + t27; t27;];
taug  = t1;
