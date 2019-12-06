% Calculate Gravitation load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:57
% EndTime: 2019-12-05 17:39:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (93->32), mult. (106->38), div. (0->0), fcn. (80->8), ass. (0->14)
t36 = sin(qJ(1));
t37 = cos(qJ(1));
t22 = g(1) * t36 - g(2) * t37;
t33 = pkin(8) + qJ(4);
t28 = qJ(5) + t33;
t24 = sin(t28);
t25 = cos(t28);
t39 = (g(3) * t24 - t22 * t25) * MDP(23) + (g(3) * t25 + t22 * t24) * MDP(24);
t38 = t37 * pkin(1) + t36 * qJ(2);
t23 = g(1) * t37 + g(2) * t36;
t30 = t37 * qJ(2);
t27 = cos(t33);
t26 = sin(t33);
t1 = [(-g(1) * (-t36 * pkin(1) + t30) - g(2) * t38) * MDP(6) + (-g(1) * (t30 + (-pkin(1) - qJ(3)) * t36) - g(2) * (t37 * qJ(3) + t38)) * MDP(10) + (MDP(2) - MDP(4) + MDP(9)) * t22 + (-t26 * MDP(16) - t27 * MDP(17) - MDP(23) * t24 - MDP(24) * t25 - MDP(7) * sin(pkin(8)) - MDP(8) * cos(pkin(8)) + MDP(3) - MDP(5)) * t23; (-MDP(10) - MDP(6)) * t22; -t23 * MDP(10); (g(3) * t26 - t22 * t27) * MDP(16) + (g(3) * t27 + t22 * t26) * MDP(17) + t39; t39;];
taug = t1;
