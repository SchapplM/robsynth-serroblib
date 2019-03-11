% Calculate Gravitation load on the joints for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:27
% EndTime: 2019-03-09 03:42:28
% DurationCPUTime: 0.30s
% Computational Cost: add. (286->68), mult. (257->105), div. (0->0), fcn. (242->12), ass. (0->37)
t66 = qJ(1) + pkin(10);
t61 = sin(t66);
t63 = cos(t66);
t80 = g(1) * t63 + g(2) * t61;
t90 = MDP(11) - MDP(14);
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t55 = -g(3) * t71 + t69 * t80;
t87 = g(3) * t69;
t85 = t61 * t71;
t84 = t63 * t71;
t67 = sin(pkin(11));
t83 = t67 * t71;
t68 = cos(pkin(11));
t82 = t68 * t71;
t65 = pkin(11) + qJ(5);
t64 = qJ(6) + t65;
t58 = sin(t64);
t59 = cos(t64);
t47 = t58 * t85 + t63 * t59;
t48 = t63 * t58 - t59 * t85;
t49 = -t58 * t84 + t61 * t59;
t50 = t61 * t58 + t59 * t84;
t81 = (-g(1) * t49 + g(2) * t47 + t58 * t87) * MDP(28) + (g(1) * t50 - g(2) * t48 + t59 * t87) * MDP(29);
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t78 = g(1) * t70 - g(2) * t72;
t77 = t71 * pkin(3) + t69 * qJ(4);
t75 = pkin(2) + t77;
t74 = t78 * pkin(1);
t62 = cos(t65);
t60 = sin(t65);
t54 = t61 * t60 + t62 * t84;
t53 = -t60 * t84 + t61 * t62;
t52 = t63 * t60 - t62 * t85;
t51 = t60 * t85 + t63 * t62;
t1 = [t78 * MDP(2) + (g(1) * t72 + g(2) * t70) * MDP(3) + MDP(4) * t74 + (-g(1) * (-t61 * t82 + t63 * t67) - g(2) * (t61 * t67 + t63 * t82)) * MDP(12) + (-g(1) * (t61 * t83 + t63 * t68) - g(2) * (t61 * t68 - t63 * t83)) * MDP(13) + (t74 + (-g(1) * pkin(7) - g(2) * t75) * t63 + (-g(2) * pkin(7) + g(1) * t75) * t61) * MDP(15) + (-g(1) * t52 - g(2) * t54) * MDP(21) + (-g(1) * t51 - g(2) * t53) * MDP(22) + (-g(1) * t48 - g(2) * t50) * MDP(28) + (-g(1) * t47 - g(2) * t49) * MDP(29) + (t71 * MDP(10) - t90 * t69) * (g(1) * t61 - g(2) * t63); (-MDP(15) - MDP(4)) * g(3); (-g(3) * t77 + t80 * (pkin(3) * t69 - qJ(4) * t71)) * MDP(15) + t90 * (t71 * t80 + t87) + (MDP(12) * t68 - MDP(13) * t67 + MDP(21) * t62 - MDP(22) * t60 + MDP(28) * t59 - MDP(29) * t58 + MDP(10)) * t55; -t55 * MDP(15); (-g(1) * t53 + g(2) * t51 + t60 * t87) * MDP(21) + (g(1) * t54 - g(2) * t52 + t62 * t87) * MDP(22) + t81; t81;];
taug  = t1;
