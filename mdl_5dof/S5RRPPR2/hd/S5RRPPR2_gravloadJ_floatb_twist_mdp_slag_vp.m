% Calculate Gravitation load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:36
% EndTime: 2020-01-03 11:57:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (177->43), mult. (139->66), div. (0->0), fcn. (118->10), ass. (0->29)
t61 = sin(pkin(9));
t74 = g(1) * t61;
t62 = cos(pkin(9));
t63 = sin(qJ(5));
t73 = t62 * t63;
t65 = cos(qJ(5));
t72 = t62 * t65;
t60 = qJ(1) + qJ(2);
t55 = pkin(8) + t60;
t51 = sin(t55);
t52 = cos(t55);
t57 = cos(t60);
t54 = pkin(2) * t57;
t71 = t52 * pkin(3) + t51 * qJ(4) + t54;
t56 = sin(t60);
t53 = pkin(2) * t56;
t70 = t51 * pkin(3) - t52 * qJ(4) + t53;
t69 = g(2) * t52 + g(3) * t51;
t68 = -g(2) * t57 - g(3) * t56;
t41 = -t51 * t73 - t52 * t65;
t42 = t51 * t72 - t52 * t63;
t43 = -t51 * t65 + t52 * t73;
t44 = t51 * t63 + t52 * t72;
t67 = (g(2) * t43 - g(3) * t41) * MDP(18) + (-g(2) * t44 - g(3) * t42) * MDP(17) + (-g(2) * t51 + g(3) * t52) * MDP(10) + (g(2) * t56 - g(3) * t57) * MDP(6) + t68 * MDP(5) + (-t62 * MDP(8) + t61 * MDP(9)) * t69;
t66 = cos(qJ(1));
t64 = sin(qJ(1));
t59 = t66 * pkin(1);
t58 = t64 * pkin(1);
t1 = [(-g(2) * t66 - g(3) * t64) * MDP(2) + (g(2) * t64 - g(3) * t66) * MDP(3) + (-g(2) * (t54 + t59) - g(3) * (t53 + t58)) * MDP(7) + (-g(2) * (t59 + t71) - g(3) * (t58 + t70)) * MDP(11) + t67; (-g(2) * t71 - g(3) * t70) * MDP(11) + t68 * MDP(7) * pkin(2) + t67; (-MDP(11) - MDP(7)) * g(1); t69 * MDP(11); (-g(2) * t41 - g(3) * t43 + t63 * t74) * MDP(17) + (g(2) * t42 - g(3) * t44 + t65 * t74) * MDP(18);];
taug = t1;
