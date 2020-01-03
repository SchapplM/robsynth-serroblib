% Calculate Gravitation load on the joints for
% S4RRPR6
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:54
% DurationCPUTime: 0.09s
% Computational Cost: add. (73->24), mult. (86->32), div. (0->0), fcn. (65->6), ass. (0->13)
t22 = qJ(2) + pkin(7) + qJ(4);
t19 = sin(t22);
t20 = cos(t22);
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t29 = g(1) * t27 + g(2) * t25;
t30 = (-g(3) * t20 + t29 * t19) * MDP(18) + (g(3) * t19 + t29 * t20) * MDP(19);
t17 = g(1) * t25 - g(2) * t27;
t26 = cos(qJ(2));
t24 = sin(qJ(2));
t23 = -qJ(3) - pkin(5);
t21 = t26 * pkin(2) + pkin(1);
t1 = [(-g(1) * (-t25 * t21 - t27 * t23) - g(2) * (t27 * t21 - t25 * t23)) * MDP(12) + (MDP(3) - MDP(11)) * t29 + (-MDP(10) * t24 + MDP(18) * t20 - MDP(19) * t19 + MDP(9) * t26 + MDP(2)) * t17; (g(3) * t24 + t29 * t26) * MDP(10) + t30 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t26 + t29 * t24); -t17 * MDP(12); t30;];
taug = t1;
