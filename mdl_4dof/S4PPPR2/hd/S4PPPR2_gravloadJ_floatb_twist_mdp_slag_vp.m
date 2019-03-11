% Calculate Gravitation load on the joints for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:18
% EndTime: 2019-03-08 18:10:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->8), mult. (22->15), div. (0->0), fcn. (18->4), ass. (0->8)
t13 = -MDP(2) - MDP(3);
t12 = cos(qJ(4));
t11 = sin(qJ(4));
t10 = cos(pkin(5));
t9 = sin(pkin(5));
t8 = t10 * t11 - t9 * t12;
t7 = -t10 * t12 - t9 * t11;
t1 = [(-MDP(1) + t13) * g(2); t13 * g(3); (-g(1) * t9 + g(2) * t10) * MDP(3); (g(1) * t8 - g(2) * t7) * MDP(5) + (-g(1) * t7 - g(2) * t8) * MDP(6);];
taug  = t1;
