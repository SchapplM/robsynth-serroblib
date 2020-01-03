% Calculate Gravitation load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:28
% EndTime: 2019-12-31 19:58:29
% DurationCPUTime: 0.32s
% Computational Cost: add. (133->54), mult. (191->75), div. (0->0), fcn. (163->8), ass. (0->30)
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t46 = g(1) * t62 + g(2) * t59;
t54 = qJ(2) + pkin(8);
t50 = sin(t54);
t51 = cos(t54);
t79 = -g(3) * t51 + t46 * t50;
t74 = g(3) * t50;
t57 = sin(qJ(4));
t72 = t59 * t57;
t60 = cos(qJ(4));
t71 = t59 * t60;
t70 = t62 * t57;
t69 = t62 * t60;
t56 = -qJ(3) - pkin(6);
t67 = pkin(4) * t57 - t56;
t45 = g(1) * t59 - g(2) * t62;
t48 = t60 * pkin(4) + pkin(3);
t55 = -qJ(5) - pkin(7);
t66 = t51 * t48 - t50 * t55;
t43 = -t51 * t70 + t71;
t41 = t51 * t72 + t69;
t61 = cos(qJ(2));
t58 = sin(qJ(2));
t52 = t61 * pkin(2);
t49 = t52 + pkin(1);
t47 = t62 * t49;
t44 = t51 * t69 + t72;
t42 = -t51 * t71 + t70;
t1 = [(-g(1) * (-t59 * t49 - t62 * t56) - g(2) * (-t59 * t56 + t47)) * MDP(12) + (-g(1) * t42 - g(2) * t44) * MDP(18) + (-g(1) * t41 - g(2) * t43) * MDP(19) + (-g(2) * t47 + (-g(1) * t67 - g(2) * t66) * t62 + (-g(1) * (-t49 - t66) - g(2) * t67) * t59) * MDP(21) + (MDP(3) - MDP(11)) * t46 + (-t58 * MDP(10) + t50 * MDP(20) + t61 * MDP(9) + MDP(2)) * t45; (g(3) * t58 + t46 * t61) * MDP(10) + (-t46 * t51 - t74) * MDP(20) + (-g(3) * (t52 + t66) + t46 * (pkin(2) * t58 + t48 * t50 + t51 * t55)) * MDP(21) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t61 + t46 * t58) + (MDP(18) * t60 - MDP(19) * t57) * t79; (-MDP(12) - MDP(21)) * t45; (g(1) * t44 - g(2) * t42 + t60 * t74) * MDP(19) + (pkin(4) * MDP(21) + MDP(18)) * (-g(1) * t43 + g(2) * t41 + t57 * t74); -t79 * MDP(21);];
taug = t1;
