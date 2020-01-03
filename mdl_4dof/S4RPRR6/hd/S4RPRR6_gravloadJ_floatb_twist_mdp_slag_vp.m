% Calculate Gravitation load on the joints for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (80->24), mult. (85->32), div. (0->0), fcn. (66->8), ass. (0->12)
t24 = pkin(7) + qJ(3);
t23 = qJ(4) + t24;
t19 = sin(t23);
t20 = cos(t23);
t27 = sin(qJ(1));
t28 = cos(qJ(1));
t29 = g(1) * t28 + g(2) * t27;
t30 = (-g(3) * t20 + t29 * t19) * MDP(20) + (g(3) * t19 + t29 * t20) * MDP(21);
t17 = g(1) * t27 - g(2) * t28;
t22 = cos(t24);
t21 = sin(t24);
t1 = [(-g(1) * (-t27 * pkin(1) + t28 * qJ(2)) - g(2) * (t28 * pkin(1) + t27 * qJ(2))) * MDP(7) + (MDP(3) - MDP(6)) * t29 + (t22 * MDP(13) - t21 * MDP(14) + MDP(20) * t20 - MDP(21) * t19 + MDP(4) * cos(pkin(7)) - MDP(5) * sin(pkin(7)) + MDP(2)) * t17; -t17 * MDP(7); (-g(3) * t22 + t29 * t21) * MDP(13) + (g(3) * t21 + t29 * t22) * MDP(14) + t30; t30;];
taug = t1;
