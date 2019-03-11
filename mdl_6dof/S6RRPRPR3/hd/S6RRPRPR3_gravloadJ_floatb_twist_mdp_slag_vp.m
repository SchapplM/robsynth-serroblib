% Calculate Gravitation load on the joints for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:13
% EndTime: 2019-03-09 10:19:14
% DurationCPUTime: 0.38s
% Computational Cost: add. (256->71), mult. (271->99), div. (0->0), fcn. (247->10), ass. (0->42)
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t65 = g(1) * t84 + g(2) * t81;
t76 = qJ(2) + pkin(10);
t71 = sin(t76);
t72 = cos(t76);
t87 = -g(3) * t72 + t65 * t71;
t102 = g(3) * t71;
t73 = qJ(4) + pkin(11) + qJ(6);
t67 = sin(t73);
t100 = t81 * t67;
t68 = cos(t73);
t99 = t81 * t68;
t79 = sin(qJ(4));
t98 = t81 * t79;
t82 = cos(qJ(4));
t97 = t81 * t82;
t96 = t84 * t67;
t95 = t84 * t68;
t94 = t84 * t79;
t93 = t84 * t82;
t56 = t72 * t100 + t95;
t57 = -t72 * t99 + t96;
t58 = -t72 * t96 + t99;
t59 = t72 * t95 + t100;
t92 = (-g(1) * t58 + g(2) * t56 + t67 * t102) * MDP(27) + (g(1) * t59 - g(2) * t57 + t68 * t102) * MDP(28);
t78 = -qJ(3) - pkin(7);
t90 = pkin(4) * t79 - t78;
t64 = g(1) * t81 - g(2) * t84;
t69 = t82 * pkin(4) + pkin(3);
t77 = -qJ(5) - pkin(8);
t89 = t72 * t69 - t71 * t77;
t62 = -t72 * t94 + t97;
t60 = t72 * t98 + t93;
t83 = cos(qJ(2));
t80 = sin(qJ(2));
t74 = t83 * pkin(2);
t70 = t74 + pkin(1);
t66 = t84 * t70;
t63 = t72 * t93 + t98;
t61 = -t72 * t97 + t94;
t1 = [(-g(1) * (-t81 * t70 - t84 * t78) - g(2) * (-t81 * t78 + t66)) * MDP(12) + (-g(1) * t61 - g(2) * t63) * MDP(18) + (-g(1) * t60 - g(2) * t62) * MDP(19) + (-g(2) * t66 + (-g(1) * t90 - g(2) * t89) * t84 + (-g(1) * (-t70 - t89) - g(2) * t90) * t81) * MDP(21) + (-g(1) * t57 - g(2) * t59) * MDP(27) + (-g(1) * t56 - g(2) * t58) * MDP(28) + (MDP(3) - MDP(11)) * t65 + (-t80 * MDP(10) + t71 * MDP(20) + t83 * MDP(9) + MDP(2)) * t64; (g(3) * t80 + t65 * t83) * MDP(10) + (-t65 * t72 - t102) * MDP(20) + (-g(3) * (t74 + t89) + t65 * (pkin(2) * t80 + t69 * t71 + t72 * t77)) * MDP(21) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t83 + t65 * t80) + (MDP(18) * t82 - MDP(19) * t79 + MDP(27) * t68 - MDP(28) * t67) * t87; (-MDP(12) - MDP(21)) * t64; (g(1) * t63 - g(2) * t61 + t82 * t102) * MDP(19) + t92 + (pkin(4) * MDP(21) + MDP(18)) * (-g(1) * t62 + g(2) * t60 + t79 * t102); -t87 * MDP(21); t92;];
taug  = t1;
