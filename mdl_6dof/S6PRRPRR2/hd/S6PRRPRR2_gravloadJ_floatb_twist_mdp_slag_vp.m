% Calculate Gravitation load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:22
% EndTime: 2019-03-08 22:00:24
% DurationCPUTime: 0.54s
% Computational Cost: add. (330->94), mult. (552->168), div. (0->0), fcn. (656->14), ass. (0->44)
t78 = sin(pkin(6));
t110 = g(3) * t78;
t75 = qJ(3) + pkin(12);
t72 = cos(t75);
t76 = qJ(5) + qJ(6);
t73 = sin(t76);
t109 = t72 * t73;
t74 = cos(t76);
t108 = t72 * t74;
t81 = sin(qJ(5));
t107 = t72 * t81;
t84 = cos(qJ(5));
t106 = t72 * t84;
t86 = cos(qJ(2));
t105 = t72 * t86;
t77 = sin(pkin(11));
t104 = t77 * t78;
t79 = cos(pkin(11));
t103 = t78 * t79;
t82 = sin(qJ(3));
t102 = t78 * t82;
t83 = sin(qJ(2));
t101 = t78 * t83;
t85 = cos(qJ(3));
t100 = t78 * t85;
t99 = t78 * t86;
t98 = t81 * t86;
t97 = t84 * t86;
t95 = cos(pkin(6));
t94 = t83 * t95;
t65 = t77 * t86 + t79 * t94;
t71 = sin(t75);
t59 = -t103 * t71 + t65 * t72;
t67 = -t77 * t94 + t79 * t86;
t61 = t104 * t71 + t67 * t72;
t63 = t101 * t72 + t71 * t95;
t93 = t86 * t95;
t64 = t77 * t83 - t79 * t93;
t66 = t77 * t93 + t79 * t83;
t96 = (-g(1) * (-t61 * t73 + t66 * t74) - g(2) * (-t59 * t73 + t64 * t74) - g(3) * (-t63 * t73 - t74 * t99)) * MDP(26) + (-g(1) * (-t61 * t74 - t66 * t73) - g(2) * (-t59 * t74 - t64 * t73) - g(3) * (-t63 * t74 + t73 * t99)) * MDP(27);
t89 = -g(1) * t66 - g(2) * t64 + g(3) * t99;
t80 = -qJ(4) - pkin(8);
t70 = pkin(3) * t85 + pkin(2);
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(1) * (-t66 * t70 - t67 * t80) - g(2) * (-t64 * t70 - t65 * t80) - (t70 * t86 - t80 * t83) * t110) * MDP(13) + (-g(1) * (-t106 * t66 + t67 * t81) - g(2) * (-t106 * t64 + t65 * t81) - (t72 * t97 + t81 * t83) * t110) * MDP(19) + (-g(1) * (t107 * t66 + t67 * t84) - g(2) * (t107 * t64 + t65 * t84) - (-t72 * t98 + t83 * t84) * t110) * MDP(20) + (-g(1) * (-t108 * t66 + t67 * t73) - g(2) * (-t108 * t64 + t65 * t73) - (t105 * t74 + t73 * t83) * t110) * MDP(26) + (-g(1) * (t109 * t66 + t67 * t74) - g(2) * (t109 * t64 + t65 * t74) - (-t105 * t73 + t74 * t83) * t110) * MDP(27) + (MDP(4) - MDP(12)) * (g(1) * t67 + g(2) * t65 + g(3) * t101) + (-t85 * MDP(10) + t82 * MDP(11) - MDP(3)) * t89; (-g(1) * (-t102 * t77 - t67 * t85) - g(2) * (t102 * t79 - t65 * t85) - g(3) * (-t100 * t83 - t82 * t95)) * MDP(11) + (-MDP(19) * t84 + MDP(20) * t81 - MDP(26) * t74 + MDP(27) * t73) * (g(1) * (t104 * t72 - t67 * t71) + g(2) * (-t103 * t72 - t65 * t71) + g(3) * (-t101 * t71 + t72 * t95)) + (MDP(13) * pkin(3) + MDP(10)) * (-g(1) * (t100 * t77 - t67 * t82) - g(2) * (-t100 * t79 - t65 * t82) - g(3) * (-t101 * t82 + t85 * t95)); t89 * MDP(13); (-g(1) * (-t61 * t81 + t66 * t84) - g(2) * (-t59 * t81 + t64 * t84) - g(3) * (-t63 * t81 - t78 * t97)) * MDP(19) + (-g(1) * (-t61 * t84 - t66 * t81) - g(2) * (-t59 * t84 - t64 * t81) - g(3) * (-t63 * t84 + t78 * t98)) * MDP(20) + t96; t96;];
taug  = t1;
