% Calculate Gravitation load on the joints for
% S4PRRP5
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
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (43->25), mult. (99->39), div. (0->0), fcn. (85->6), ass. (0->14)
t23 = sin(pkin(6));
t24 = cos(pkin(6));
t32 = g(1) * t24 + g(2) * t23;
t27 = sin(qJ(2));
t29 = cos(qJ(2));
t20 = -g(3) * t29 + t32 * t27;
t36 = g(3) * t27;
t26 = sin(qJ(3));
t34 = t26 * t29;
t28 = cos(qJ(3));
t33 = t28 * t29;
t25 = -qJ(4) - pkin(5);
t22 = t28 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(3) * (t29 * t22 - t27 * t25) + t32 * (t22 * t27 + t25 * t29)) * MDP(13) + (MDP(4) - MDP(12)) * (t32 * t29 + t36) + (MDP(10) * t28 - MDP(11) * t26 + MDP(3)) * t20; (-g(1) * (-t23 * t26 - t24 * t33) - g(2) * (-t23 * t33 + t24 * t26) + t28 * t36) * MDP(11) + (pkin(3) * MDP(13) + MDP(10)) * (-g(1) * (t23 * t28 - t24 * t34) - g(2) * (-t23 * t34 - t24 * t28) + t26 * t36); -t20 * MDP(13);];
taug = t1;
