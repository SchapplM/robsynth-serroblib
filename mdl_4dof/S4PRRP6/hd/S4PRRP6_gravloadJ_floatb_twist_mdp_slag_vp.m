% Calculate Gravitation load on the joints for
% S4PRRP6
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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (66->32), mult. (165->49), div. (0->0), fcn. (157->6), ass. (0->19)
t56 = MDP(10) + MDP(12);
t55 = MDP(11) - MDP(14);
t38 = sin(pkin(6));
t39 = cos(pkin(6));
t54 = -g(1) * t39 - g(2) * t38;
t41 = sin(qJ(2));
t51 = g(3) * t41;
t40 = sin(qJ(3));
t43 = cos(qJ(2));
t50 = t40 * t43;
t42 = cos(qJ(3));
t49 = t42 * t43;
t47 = pkin(3) * t42 + qJ(4) * t40 + pkin(2);
t32 = t38 * t50 + t39 * t42;
t34 = -t38 * t42 + t39 * t50;
t28 = g(1) * t34 + g(2) * t32 + t40 * t51;
t35 = t38 * t40 + t39 * t49;
t33 = t38 * t49 - t39 * t40;
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(3) * (t41 * pkin(5) + t47 * t43) + t54 * (pkin(5) * t43 - t47 * t41)) * MDP(15) + (MDP(4) - MDP(13)) * (-t54 * t43 + t51) + (-t55 * t40 + t56 * t42 + MDP(3)) * (-g(3) * t43 - t54 * t41); (-g(1) * (-t34 * pkin(3) + t35 * qJ(4)) - g(2) * (-t32 * pkin(3) + t33 * qJ(4)) - (-pkin(3) * t40 + qJ(4) * t42) * t51) * MDP(15) + t55 * (g(1) * t35 + g(2) * t33 + t42 * t51) + t56 * t28; -t28 * MDP(15);];
taug = t1;
