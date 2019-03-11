% Calculate Gravitation load on the joints for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:00
% EndTime: 2019-03-08 18:19:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->18), mult. (34->25), div. (0->0), fcn. (18->4), ass. (0->8)
t15 = -MDP(5) - MDP(8);
t13 = cos(qJ(2));
t12 = sin(qJ(2));
t11 = qJ(2) + pkin(5);
t10 = cos(t11);
t9 = sin(t11);
t8 = g(1) * t9 - g(2) * t10;
t1 = [(-MDP(1) + t15) * g(2); (g(1) * t13 + g(2) * t12) * MDP(4) + t8 * MDP(6) + (-g(1) * t10 - g(2) * t9) * MDP(7) + (-g(1) * (-t12 * pkin(2) - t9 * pkin(3) + t10 * qJ(4)) - g(2) * (t13 * pkin(2) + t10 * pkin(3) + t9 * qJ(4))) * MDP(8) + (pkin(2) * MDP(5) + MDP(3)) * (g(1) * t12 - g(2) * t13); t15 * g(3); -t8 * MDP(8);];
taug  = t1;
