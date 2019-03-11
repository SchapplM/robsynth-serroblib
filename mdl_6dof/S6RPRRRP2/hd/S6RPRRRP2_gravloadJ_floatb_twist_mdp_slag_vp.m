% Calculate Gravitation load on the joints for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:16
% EndTime: 2019-03-09 06:01:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (268->68), mult. (274->104), div. (0->0), fcn. (251->10), ass. (0->41)
t76 = qJ(1) + pkin(10);
t70 = sin(t76);
t71 = cos(t76);
t91 = g(1) * t71 + g(2) * t70;
t107 = MDP(11) - MDP(26);
t79 = sin(qJ(3));
t101 = g(3) * t79;
t77 = qJ(4) + qJ(5);
t73 = cos(t77);
t72 = sin(t77);
t82 = cos(qJ(3));
t97 = t72 * t82;
t55 = t70 * t97 + t71 * t73;
t57 = t70 * t73 - t71 * t97;
t106 = -g(1) * t57 + g(2) * t55 + t72 * t101;
t59 = -g(3) * t82 + t91 * t79;
t78 = sin(qJ(4));
t67 = t78 * pkin(4) + pkin(5) * t72;
t99 = pkin(7) + t67;
t98 = t67 * t82;
t96 = t73 * t82;
t95 = t78 * t82;
t81 = cos(qJ(4));
t94 = t81 * t82;
t56 = -t70 * t96 + t71 * t72;
t58 = t70 * t72 + t71 * t96;
t93 = t106 * MDP(24) + (g(1) * t58 - g(2) * t56 + t73 * t101) * MDP(25);
t68 = t81 * pkin(4) + pkin(5) * t73;
t80 = sin(qJ(1));
t83 = cos(qJ(1));
t89 = g(1) * t80 - g(2) * t83;
t66 = pkin(3) + t68;
t75 = -qJ(6) - pkin(9) - pkin(8);
t88 = t82 * t66 - t79 * t75;
t86 = pkin(2) + t88;
t85 = t89 * pkin(1);
t64 = t70 * t78 + t71 * t94;
t63 = t70 * t81 - t71 * t95;
t62 = -t70 * t94 + t71 * t78;
t61 = t70 * t95 + t71 * t81;
t1 = [t89 * MDP(2) + (g(1) * t83 + g(2) * t80) * MDP(3) + MDP(4) * t85 + (-g(1) * t62 - g(2) * t64) * MDP(17) + (-g(1) * t61 - g(2) * t63) * MDP(18) + (-g(1) * t56 - g(2) * t58) * MDP(24) + (-g(1) * t55 - g(2) * t57) * MDP(25) + (t85 + (-g(1) * t99 - g(2) * t86) * t71 + (g(1) * t86 - g(2) * t99) * t70) * MDP(27) + (t82 * MDP(10) - t107 * t79) * (g(1) * t70 - g(2) * t71); (-MDP(27) - MDP(4)) * g(3); (-g(3) * t88 + t91 * (t66 * t79 + t75 * t82)) * MDP(27) + t107 * (t91 * t82 + t101) + (MDP(17) * t81 - MDP(18) * t78 + MDP(24) * t73 - MDP(25) * t72 + MDP(10)) * t59; (-g(1) * t63 + g(2) * t61 + t78 * t101) * MDP(17) + (g(1) * t64 - g(2) * t62 + t81 * t101) * MDP(18) + (-g(1) * (t70 * t68 - t71 * t98) - g(2) * (-t71 * t68 - t70 * t98) + t67 * t101) * MDP(27) + t93; t106 * MDP(27) * pkin(5) + t93; -t59 * MDP(27);];
taug  = t1;
