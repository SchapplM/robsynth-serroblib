% Calculate Gravitation load on the joints for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x6] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:45
% EndTime: 2019-01-31 13:16:46
% DurationCPUTime: 0.04s
% Computational Cost: add. (76->19), mult. (49->28), div. (0->0), fcn. (30->6), ass. (0->12)
t22 = qJ(1) + qJ(2);
t19 = pkin(7) + qJ(4) + t22;
t17 = sin(t19);
t18 = cos(t19);
t27 = (g(1) * t17 - g(2) * t18) * MDP(9) + (g(1) * t18 + g(2) * t17) * MDP(10);
t20 = sin(t22);
t21 = cos(t22);
t25 = g(1) * t20 - g(2) * t21;
t26 = t25 * MDP(5) + (g(1) * t21 + g(2) * t20) * MDP(6) + t27;
t24 = cos(qJ(1));
t23 = sin(qJ(1));
t1 = [(g(1) * t23 - g(2) * t24) * MDP(2) + (g(1) * t24 + g(2) * t23) * MDP(3) + (-g(1) * (-t23 * pkin(1) - pkin(2) * t20) - g(2) * (t24 * pkin(1) + pkin(2) * t21)) * MDP(7) + t26; t25 * MDP(7) * pkin(2) + t26; -g(3) * MDP(7); t27;];
taug  = t1;
