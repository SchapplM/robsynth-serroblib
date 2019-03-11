% Calculate Gravitation load on the joints for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S2RR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S2RR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:03
% EndTime: 2019-03-08 18:00:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (12->8), mult. (28->14), div. (0->0), fcn. (22->4), ass. (0->6)
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t12 = -g(1) * t8 - g(3) * t10;
t9 = cos(qJ(2));
t7 = sin(qJ(2));
t1 = [t12 * MDP(3) + (-t7 * MDP(10) + t9 * MDP(9) + MDP(2)) * (g(1) * t10 - g(3) * t8); (g(2) * t9 + t12 * t7) * MDP(9) + (-g(2) * t7 + t12 * t9) * MDP(10);];
taug  = t1;
