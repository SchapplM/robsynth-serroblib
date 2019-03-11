% Calculate Gravitation load on the joints for
% S6PRRPRR1
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
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:07
% EndTime: 2019-03-08 21:54:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (330->73), mult. (474->128), div. (0->0), fcn. (548->12), ass. (0->37)
t72 = sin(pkin(6));
t100 = g(3) * t72;
t70 = qJ(3) + pkin(12) + qJ(5);
t68 = cos(t70);
t75 = sin(qJ(6));
t99 = t68 * t75;
t78 = cos(qJ(6));
t98 = t68 * t78;
t71 = sin(pkin(11));
t97 = t71 * t72;
t73 = cos(pkin(11));
t96 = t72 * t73;
t76 = sin(qJ(3));
t95 = t72 * t76;
t77 = sin(qJ(2));
t94 = t72 * t77;
t79 = cos(qJ(3));
t93 = t72 * t79;
t80 = cos(qJ(2));
t92 = t75 * t80;
t91 = t78 * t80;
t90 = cos(pkin(6));
t89 = t77 * t90;
t88 = t80 * t90;
t62 = t71 * t80 + t73 * t89;
t67 = sin(t70);
t56 = t62 * t68 - t67 * t96;
t64 = -t71 * t89 + t73 * t80;
t58 = t64 * t68 + t67 * t97;
t60 = t90 * t67 + t68 * t94;
t87 = (g(1) * t58 + g(2) * t56 + g(3) * t60) * MDP(20) + (-t78 * MDP(26) + t75 * MDP(27) - MDP(19)) * (g(1) * (-t64 * t67 + t68 * t97) + g(2) * (-t62 * t67 - t68 * t96) + g(3) * (-t67 * t94 + t90 * t68));
t61 = t71 * t77 - t73 * t88;
t63 = t71 * t88 + t73 * t77;
t83 = -g(1) * t63 - g(2) * t61 + t80 * t100;
t74 = -qJ(4) - pkin(8);
t69 = t79 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(1) * (-t63 * t69 - t64 * t74) - g(2) * (-t61 * t69 - t62 * t74) - (t69 * t80 - t74 * t77) * t100) * MDP(13) + (-g(1) * (-t63 * t98 + t64 * t75) - g(2) * (-t61 * t98 + t62 * t75) - (t68 * t91 + t75 * t77) * t100) * MDP(26) + (-g(1) * (t63 * t99 + t64 * t78) - g(2) * (t61 * t99 + t62 * t78) - (-t68 * t92 + t77 * t78) * t100) * MDP(27) + (MDP(4) - MDP(12)) * (g(1) * t64 + g(2) * t62 + g(3) * t94) + (-t79 * MDP(10) + t76 * MDP(11) - t68 * MDP(19) + MDP(20) * t67 - MDP(3)) * t83; (-g(1) * (-t64 * t79 - t71 * t95) - g(2) * (-t62 * t79 + t73 * t95) - g(3) * (-t90 * t76 - t77 * t93)) * MDP(11) + t87 + (pkin(3) * MDP(13) + MDP(10)) * (-g(1) * (-t64 * t76 + t71 * t93) - g(2) * (-t62 * t76 - t73 * t93) - g(3) * (-t76 * t94 + t90 * t79)); t83 * MDP(13); t87; (-g(1) * (-t58 * t75 + t63 * t78) - g(2) * (-t56 * t75 + t61 * t78) - g(3) * (-t60 * t75 - t72 * t91)) * MDP(26) + (-g(1) * (-t58 * t78 - t63 * t75) - g(2) * (-t56 * t78 - t61 * t75) - g(3) * (-t60 * t78 + t72 * t92)) * MDP(27);];
taug  = t1;
