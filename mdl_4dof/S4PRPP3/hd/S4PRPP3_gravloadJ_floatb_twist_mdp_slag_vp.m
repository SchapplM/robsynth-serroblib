% Calculate Gravitation load on the joints for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:54
% EndTime: 2019-03-08 18:19:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (26->18), mult. (44->21), div. (0->0), fcn. (26->2), ass. (0->7)
t16 = sin(qJ(2));
t17 = cos(qJ(2));
t19 = t17 * pkin(2) + t16 * qJ(3);
t18 = -MDP(10) - MDP(7);
t13 = t17 * qJ(3);
t10 = g(1) * t16 - g(2) * t17;
t1 = [(-MDP(1) + t18) * g(2); (-g(1) * (-t16 * pkin(2) + t13) - g(2) * t19) * MDP(7) + (-g(1) * (t13 + (-pkin(2) - pkin(3)) * t16) - g(2) * (t17 * pkin(3) + t19)) * MDP(10) + (MDP(4) - MDP(6) - MDP(9)) * (g(1) * t17 + g(2) * t16) + (MDP(3) + MDP(5) + MDP(8)) * t10; t18 * t10; g(3) * MDP(10);];
taug  = t1;
