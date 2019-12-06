% Calculate Gravitation load on the joints for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:36
% EndTime: 2019-12-05 15:51:37
% DurationCPUTime: 0.31s
% Computational Cost: add. (173->59), mult. (456->114), div. (0->0), fcn. (569->12), ass. (0->33)
t69 = sin(pkin(10));
t72 = cos(pkin(10));
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t65 = t77 * t69 - t80 * t72;
t71 = sin(pkin(5));
t96 = g(3) * t71;
t76 = sin(qJ(4));
t95 = t71 * t76;
t79 = cos(qJ(4));
t94 = t71 * t79;
t74 = cos(pkin(5));
t93 = t74 * t77;
t92 = t74 * t80;
t75 = sin(qJ(5));
t91 = t75 * t79;
t78 = cos(qJ(5));
t89 = t78 * t79;
t85 = t80 * t69 + t77 * t72;
t64 = t85 * t74;
t70 = sin(pkin(9));
t73 = cos(pkin(9));
t87 = t73 * t64 - t70 * t65;
t86 = -t70 * t64 - t73 * t65;
t84 = t65 * t74;
t63 = t85 * t71;
t62 = t65 * t71;
t60 = t63 * t79 + t74 * t76;
t57 = t70 * t84 - t73 * t85;
t54 = -t70 * t85 - t73 * t84;
t52 = t70 * t95 + t79 * t86;
t50 = -t73 * t95 + t79 * t87;
t1 = [(-MDP(1) - MDP(5)) * g(3); (-g(1) * (t70 * t93 - t73 * t80) - g(2) * (-t70 * t80 - t73 * t93) + t77 * t96) * MDP(4) + (-g(1) * (t57 * t89 + t75 * t86) - g(2) * (t54 * t89 + t75 * t87) - g(3) * (-t62 * t89 + t63 * t75)) * MDP(18) + (-g(1) * (-t57 * t91 + t78 * t86) - g(2) * (-t54 * t91 + t78 * t87) - g(3) * (t62 * t91 + t63 * t78)) * MDP(19) + (-t79 * MDP(11) + MDP(12) * t76) * (g(1) * t57 + g(2) * t54 - g(3) * t62) + (pkin(2) * MDP(5) + MDP(3)) * (-g(1) * (-t70 * t92 - t73 * t77) - g(2) * (-t70 * t77 + t73 * t92) - t80 * t96); (-g(3) * t74 + (-g(1) * t70 + g(2) * t73) * t71) * MDP(5); (g(1) * t52 + g(2) * t50 + g(3) * t60) * MDP(12) + (-MDP(18) * t78 + MDP(19) * t75 - MDP(11)) * (g(1) * (t70 * t94 - t76 * t86) + g(2) * (-t73 * t94 - t76 * t87) + g(3) * (-t63 * t76 + t74 * t79)); (-g(1) * (-t52 * t75 - t57 * t78) - g(2) * (-t50 * t75 - t54 * t78) - g(3) * (-t60 * t75 + t62 * t78)) * MDP(18) + (-g(1) * (-t52 * t78 + t57 * t75) - g(2) * (-t50 * t78 + t54 * t75) - g(3) * (-t60 * t78 - t62 * t75)) * MDP(19);];
taug = t1;
