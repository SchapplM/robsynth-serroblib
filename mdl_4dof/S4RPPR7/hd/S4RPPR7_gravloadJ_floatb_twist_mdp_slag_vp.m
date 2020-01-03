% Calculate Gravitation load on the joints for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:41
% EndTime: 2019-12-31 16:41:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (46->25), mult. (72->30), div. (0->0), fcn. (52->6), ass. (0->10)
t26 = sin(qJ(1));
t27 = cos(qJ(1));
t28 = t27 * pkin(1) + t26 * qJ(2);
t16 = g(1) * t27 + g(2) * t26;
t15 = g(1) * t26 - g(2) * t27;
t23 = pkin(6) + qJ(4);
t20 = t27 * qJ(2);
t18 = cos(t23);
t17 = sin(t23);
t1 = [(-g(1) * (-t26 * pkin(1) + t20) - g(2) * t28) * MDP(6) + (-g(1) * (t20 + (-pkin(1) - qJ(3)) * t26) - g(2) * (t27 * qJ(3) + t28)) * MDP(10) + (MDP(2) - MDP(4) + MDP(9)) * t15 + (-t17 * MDP(16) - t18 * MDP(17) - MDP(7) * sin(pkin(6)) - MDP(8) * cos(pkin(6)) + MDP(3) - MDP(5)) * t16; (-MDP(10) - MDP(6)) * t15; -t16 * MDP(10); (g(3) * t17 - t15 * t18) * MDP(16) + (g(3) * t18 + t15 * t17) * MDP(17);];
taug = t1;
