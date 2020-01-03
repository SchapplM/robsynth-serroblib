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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:25:49
% EndTime: 2020-01-03 11:25:51
% DurationCPUTime: 0.18s
% Computational Cost: add. (112->44), mult. (132->67), div. (0->0), fcn. (111->8), ass. (0->24)
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t58 = cos(qJ(4));
t76 = -t53 * (-qJ(5) - pkin(6)) + (t58 * pkin(4) + pkin(3)) * t54;
t74 = g(1) * t53;
t52 = qJ(1) + pkin(7);
t48 = sin(t52);
t56 = sin(qJ(4));
t70 = t48 * t56;
t68 = t54 * t56;
t67 = t54 * t58;
t57 = sin(qJ(1));
t66 = t57 * pkin(1) + t48 * pkin(2);
t65 = MDP(17) + MDP(8);
t49 = cos(t52);
t59 = cos(qJ(1));
t63 = t59 * pkin(1) + t49 * pkin(2) + t48 * qJ(3);
t62 = g(2) * t49 + g(3) * t48;
t61 = -g(2) * t48 + g(3) * t49;
t40 = -t48 * t58 + t49 * t68;
t38 = -t48 * t68 - t49 * t58;
t41 = t49 * t67 + t70;
t39 = t48 * t67 - t49 * t56;
t1 = [(g(2) * t57 - g(3) * t59) * MDP(3) + t61 * MDP(7) + (-g(2) * t63 - g(3) * (-t49 * qJ(3) + t66)) * MDP(8) + (-g(2) * t41 - g(3) * t39) * MDP(14) + (g(2) * t40 - g(3) * t38) * MDP(15) + (-g(2) * (pkin(4) * t70 + t63) - g(3) * (t76 * t48 + t66) + (-g(2) * t76 - g(3) * (-pkin(4) * t56 - qJ(3))) * t49) * MDP(17) + (MDP(4) * pkin(1) + MDP(2)) * (-g(2) * t59 - g(3) * t57) + (-t54 * MDP(5) + (MDP(6) - MDP(16)) * t53) * t62; (-MDP(4) - t65) * g(1); t65 * t62; (g(2) * t39 - g(3) * t41 + t58 * t74) * MDP(15) + (MDP(17) * pkin(4) + MDP(14)) * (-g(2) * t38 - g(3) * t40 + t56 * t74); (g(1) * t54 + t61 * t53) * MDP(17);];
taug = t1;
