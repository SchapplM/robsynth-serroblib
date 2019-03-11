% Calculate Gravitation load on the joints for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:36
% EndTime: 2019-03-09 08:52:38
% DurationCPUTime: 0.41s
% Computational Cost: add. (284->78), mult. (282->113), div. (0->0), fcn. (262->12), ass. (0->43)
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t65 = g(1) * t85 + g(2) * t83;
t78 = qJ(2) + pkin(10);
t71 = sin(t78);
t73 = cos(t78);
t88 = -g(3) * t73 + t65 * t71;
t102 = g(3) * t71;
t100 = t73 * t85;
t79 = sin(pkin(11));
t99 = t79 * t85;
t80 = cos(pkin(11));
t98 = t80 * t85;
t77 = pkin(11) + qJ(5);
t74 = qJ(6) + t77;
t67 = sin(t74);
t97 = t83 * t67;
t68 = cos(t74);
t96 = t83 * t68;
t70 = sin(t77);
t95 = t83 * t70;
t72 = cos(t77);
t94 = t83 * t72;
t93 = t83 * t79;
t92 = t83 * t80;
t56 = t68 * t85 + t73 * t97;
t57 = t67 * t85 - t73 * t96;
t58 = -t67 * t100 + t96;
t59 = t68 * t100 + t97;
t91 = (-g(1) * t58 + g(2) * t56 + t67 * t102) * MDP(29) + (g(1) * t59 - g(2) * t57 + t68 * t102) * MDP(30);
t64 = g(1) * t83 - g(2) * t85;
t90 = pkin(3) * t73 + qJ(4) * t71;
t84 = cos(qJ(2));
t82 = sin(qJ(2));
t81 = -qJ(3) - pkin(7);
t75 = t84 * pkin(2);
t69 = t75 + pkin(1);
t66 = t85 * t69;
t63 = t72 * t100 + t95;
t62 = -t70 * t100 + t94;
t61 = t70 * t85 - t73 * t94;
t60 = t72 * t85 + t73 * t95;
t1 = [(-g(1) * (-t83 * t69 - t81 * t85) - g(2) * (-t83 * t81 + t66)) * MDP(12) + (-g(1) * (-t73 * t92 + t99) - g(2) * (t73 * t98 + t93)) * MDP(13) + (-g(1) * (t73 * t93 + t98) - g(2) * (-t73 * t99 + t92)) * MDP(14) + (-g(2) * t66 + (g(1) * t81 - g(2) * t90) * t85 + (-g(1) * (-t69 - t90) + g(2) * t81) * t83) * MDP(16) + (-g(1) * t61 - g(2) * t63) * MDP(22) + (-g(1) * t60 - g(2) * t62) * MDP(23) + (-g(1) * t57 - g(2) * t59) * MDP(29) + (-g(1) * t56 - g(2) * t58) * MDP(30) + (MDP(3) - MDP(11)) * t65 + (-t82 * MDP(10) + t71 * MDP(15) + t84 * MDP(9) + MDP(2)) * t64; (g(3) * t82 + t65 * t84) * MDP(10) + (-t65 * t73 - t102) * MDP(15) + (-g(3) * (t75 + t90) + t65 * (pkin(2) * t82 + pkin(3) * t71 - qJ(4) * t73)) * MDP(16) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t84 + t65 * t82) + (t80 * MDP(13) - t79 * MDP(14) + t72 * MDP(22) - t70 * MDP(23) + t68 * MDP(29) - t67 * MDP(30)) * t88; (-MDP(12) - MDP(16)) * t64; -t88 * MDP(16); (-g(1) * t62 + g(2) * t60 + t70 * t102) * MDP(22) + (g(1) * t63 - g(2) * t61 + t72 * t102) * MDP(23) + t91; t91;];
taug  = t1;
