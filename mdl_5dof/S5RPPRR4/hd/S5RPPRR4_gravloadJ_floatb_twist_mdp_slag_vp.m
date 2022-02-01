% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:58
% EndTime: 2022-01-23 09:16:59
% DurationCPUTime: 0.16s
% Computational Cost: add. (150->54), mult. (164->84), div. (0->0), fcn. (163->10), ass. (0->40)
t65 = sin(pkin(8));
t84 = g(3) * t65;
t63 = pkin(9) + qJ(4);
t59 = qJ(5) + t63;
t55 = sin(t59);
t68 = sin(qJ(1));
t83 = t68 * t55;
t56 = cos(t59);
t82 = t68 * t56;
t57 = sin(t63);
t81 = t68 * t57;
t58 = cos(t63);
t80 = t68 * t58;
t64 = sin(pkin(9));
t79 = t68 * t64;
t66 = cos(pkin(9));
t78 = t68 * t66;
t69 = cos(qJ(1));
t77 = t69 * t55;
t76 = t69 * t56;
t75 = t69 * t57;
t74 = t69 * t58;
t73 = t69 * t64;
t72 = t69 * t66;
t67 = cos(pkin(8));
t44 = t67 * t83 + t76;
t45 = -t67 * t82 + t77;
t46 = -t67 * t77 + t82;
t47 = t67 * t76 + t83;
t71 = (-g(1) * t46 + g(2) * t44 + t55 * t84) * MDP(22) + (g(1) * t47 - g(2) * t45 + t56 * t84) * MDP(23);
t70 = g(1) * t69 + g(2) * t68;
t53 = g(1) * t68 - g(2) * t69;
t61 = t69 * qJ(2);
t60 = t68 * qJ(2);
t52 = pkin(2) * t67 + t65 * qJ(3) + pkin(1);
t51 = t67 * t74 + t81;
t50 = -t67 * t75 + t80;
t49 = -t67 * t80 + t75;
t48 = t67 * t81 + t74;
t1 = [(-g(1) * (-t68 * pkin(1) + t61) - g(2) * (t69 * pkin(1) + t60)) * MDP(6) + (-g(1) * (-t67 * t78 + t73) - g(2) * (t67 * t72 + t79)) * MDP(7) + (-g(1) * (t67 * t79 + t72) - g(2) * (-t67 * t73 + t78)) * MDP(8) + (-g(1) * (-t52 * t68 + t61) - g(2) * (t52 * t69 + t60)) * MDP(9) + (-g(1) * t49 - g(2) * t51) * MDP(15) + (-g(1) * t48 - g(2) * t50) * MDP(16) + (-g(1) * t45 - g(2) * t47) * MDP(22) + (-g(1) * t44 - g(2) * t46) * MDP(23) + (MDP(3) - MDP(5)) * t70 + (t67 * MDP(4) + MDP(2)) * t53; (-MDP(6) - MDP(9)) * t53; (g(3) * t67 - t70 * t65) * MDP(9); (-g(1) * t50 + g(2) * t48 + t57 * t84) * MDP(15) + (g(1) * t51 - g(2) * t49 + t58 * t84) * MDP(16) + t71; t71;];
taug = t1;
