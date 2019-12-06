% Calculate Gravitation load on the joints for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:13
% EndTime: 2019-12-05 16:40:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (128->28), mult. (90->36), div. (0->0), fcn. (63->6), ass. (0->16)
t33 = pkin(8) + qJ(2);
t32 = qJ(3) + t33;
t27 = sin(t32);
t28 = cos(t32);
t36 = cos(qJ(4));
t29 = t36 * pkin(4) + pkin(3);
t34 = -qJ(5) - pkin(7);
t40 = -t27 * t34 + t28 * t29;
t24 = g(1) * t27 - g(2) * t28;
t25 = g(1) * t28 + g(2) * t27;
t35 = sin(qJ(4));
t39 = (-MDP(15) + MDP(7)) * t25 + (t36 * MDP(13) - t35 * MDP(14) + MDP(6)) * t24;
t38 = -t27 * t29 - t28 * t34;
t31 = cos(t33);
t30 = sin(t33);
t1 = [(-MDP(1) - MDP(16)) * g(3); (g(1) * t30 - g(2) * t31) * MDP(3) + (g(1) * t31 + g(2) * t30) * MDP(4) + (-g(1) * (-pkin(2) * t30 + t38) - g(2) * (pkin(2) * t31 + t40)) * MDP(16) + t39; (-g(1) * t38 - g(2) * t40) * MDP(16) + t39; (g(3) * t35 + t25 * t36) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (-g(3) * t36 + t25 * t35); -t24 * MDP(16);];
taug = t1;
