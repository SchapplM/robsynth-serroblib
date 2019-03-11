% Calculate Gravitation load on the joints for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:00
% EndTime: 2019-03-08 18:23:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (88->21), mult. (51->26), div. (0->0), fcn. (32->4), ass. (0->11)
t28 = pkin(6) + qJ(2);
t27 = qJ(3) + t28;
t23 = sin(t27);
t24 = cos(t27);
t31 = t24 * pkin(3) + t23 * qJ(4);
t30 = -t23 * pkin(3) + t24 * qJ(4);
t18 = g(1) * t23 - g(2) * t24;
t29 = (MDP(7) - MDP(9)) * (g(1) * t24 + g(2) * t23) + (MDP(6) + MDP(8)) * t18;
t26 = cos(t28);
t25 = sin(t28);
t1 = [(-MDP(1) - MDP(10)) * g(3); (g(1) * t25 - g(2) * t26) * MDP(3) + (g(1) * t26 + g(2) * t25) * MDP(4) + (-g(1) * (-pkin(2) * t25 + t30) - g(2) * (pkin(2) * t26 + t31)) * MDP(10) + t29; (-g(1) * t30 - g(2) * t31) * MDP(10) + t29; -t18 * MDP(10);];
taug  = t1;
