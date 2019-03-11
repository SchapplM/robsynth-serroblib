% Calculate Gravitation load on the joints for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:08
% EndTime: 2019-03-08 18:08:08
% DurationCPUTime: 0.03s
% Computational Cost: add. (53->13), mult. (36->18), div. (0->0), fcn. (24->6), ass. (0->11)
t22 = qJ(1) + qJ(2);
t21 = qJ(3) + t22;
t17 = sin(t21);
t18 = cos(t21);
t26 = (g(1) * t17 - g(2) * t18) * MDP(8) + (g(1) * t18 + g(2) * t17) * MDP(9);
t19 = sin(t22);
t20 = cos(t22);
t25 = (g(1) * t19 - g(2) * t20) * MDP(5) + (g(1) * t20 + g(2) * t19) * MDP(6) + t26;
t24 = cos(qJ(1));
t23 = sin(qJ(1));
t1 = [(g(1) * t23 - g(2) * t24) * MDP(2) + (g(1) * t24 + g(2) * t23) * MDP(3) + t25; t25; t26;];
taug  = t1;
