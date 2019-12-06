% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:25
% EndTime: 2019-12-05 17:36:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (112->46), mult. (132->68), div. (0->0), fcn. (111->8), ass. (0->24)
t47 = sin(qJ(4));
t64 = pkin(4) * t47;
t44 = sin(pkin(8));
t63 = g(1) * t44;
t50 = cos(qJ(1));
t60 = t50 * pkin(1);
t45 = cos(pkin(8));
t59 = t45 * t47;
t49 = cos(qJ(4));
t58 = t45 * t49;
t57 = MDP(17) + MDP(8);
t43 = qJ(1) + pkin(7);
t42 = cos(t43);
t48 = sin(qJ(1));
t55 = -t48 * pkin(1) + t42 * qJ(3);
t41 = sin(t43);
t54 = g(2) * t42 + g(3) * t41;
t53 = g(2) * t41 - g(3) * t42;
t51 = (t49 * pkin(4) + pkin(3)) * t45 - t44 * (-qJ(5) - pkin(6)) + pkin(2);
t35 = -t41 * t49 + t42 * t59;
t33 = t41 * t59 + t42 * t49;
t36 = -t41 * t47 - t42 * t58;
t34 = t41 * t58 - t42 * t47;
t1 = [(-g(2) * t48 + g(3) * t50) * MDP(3) + t53 * MDP(7) + (-g(2) * (-t42 * pkin(2) - t41 * qJ(3) - t60) - g(3) * (-t41 * pkin(2) + t55)) * MDP(8) + (-g(2) * t36 + g(3) * t34) * MDP(14) + (-g(2) * t35 - g(3) * t33) * MDP(15) + (g(2) * t60 - g(3) * t55 + (g(2) * t51 - g(3) * t64) * t42 + (-g(2) * (-qJ(3) - t64) + g(3) * t51) * t41) * MDP(17) + (pkin(1) * MDP(4) + MDP(2)) * (g(2) * t50 + g(3) * t48) + (t45 * MDP(5) + (-MDP(6) + MDP(16)) * t44) * t54; (-MDP(4) - t57) * g(1); -t57 * t54; (-g(2) * t34 - g(3) * t36 + t49 * t63) * MDP(15) + (MDP(17) * pkin(4) + MDP(14)) * (-g(2) * t33 + g(3) * t35 + t47 * t63); (g(1) * t45 + t53 * t44) * MDP(17);];
taug = t1;
