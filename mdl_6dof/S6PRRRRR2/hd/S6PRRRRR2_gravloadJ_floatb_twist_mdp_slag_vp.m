% Calculate Gravitation load on the joints for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:23
% EndTime: 2019-03-09 00:45:24
% DurationCPUTime: 0.43s
% Computational Cost: add. (430->90), mult. (669->160), div. (0->0), fcn. (808->14), ass. (0->41)
t79 = sin(pkin(6));
t106 = g(3) * t79;
t76 = qJ(5) + qJ(6);
t72 = sin(t76);
t77 = qJ(3) + qJ(4);
t75 = cos(t77);
t105 = t72 * t75;
t74 = cos(t76);
t104 = t74 * t75;
t80 = sin(qJ(5));
t103 = t75 * t80;
t83 = cos(qJ(5));
t102 = t75 * t83;
t85 = cos(qJ(2));
t101 = t75 * t85;
t78 = sin(pkin(12));
t100 = t78 * t79;
t82 = sin(qJ(2));
t99 = t79 * t82;
t84 = cos(qJ(3));
t98 = t79 * t84;
t97 = t79 * t85;
t96 = t80 * t85;
t95 = t83 * t85;
t92 = cos(pkin(12));
t93 = cos(pkin(6));
t89 = t93 * t92;
t68 = t78 * t85 + t82 * t89;
t73 = sin(t77);
t90 = t79 * t92;
t62 = t68 * t75 - t73 * t90;
t91 = t78 * t93;
t70 = -t82 * t91 + t92 * t85;
t64 = t73 * t100 + t70 * t75;
t66 = t93 * t73 + t75 * t99;
t67 = t78 * t82 - t85 * t89;
t69 = t92 * t82 + t85 * t91;
t94 = (-g(1) * (-t64 * t72 + t69 * t74) - g(2) * (-t62 * t72 + t67 * t74) - g(3) * (-t66 * t72 - t74 * t97)) * MDP(31) + (-g(1) * (-t64 * t74 - t69 * t72) - g(2) * (-t62 * t74 - t67 * t72) - g(3) * (-t66 * t74 + t72 * t97)) * MDP(32);
t88 = (g(1) * t64 + g(2) * t62 + g(3) * t66) * MDP(18) + (-t83 * MDP(24) + t80 * MDP(25) - t74 * MDP(31) + t72 * MDP(32) - MDP(17)) * (g(1) * (t75 * t100 - t70 * t73) + g(2) * (-t68 * t73 - t75 * t90) + g(3) * (-t73 * t99 + t93 * t75));
t81 = sin(qJ(3));
t1 = [-g(3) * MDP(1); (g(1) * t70 + g(2) * t68 + g(3) * t99) * MDP(4) + (-g(1) * (-t69 * t102 + t70 * t80) - g(2) * (-t67 * t102 + t68 * t80) - (t75 * t95 + t80 * t82) * t106) * MDP(24) + (-g(1) * (t69 * t103 + t70 * t83) - g(2) * (t67 * t103 + t68 * t83) - (-t75 * t96 + t82 * t83) * t106) * MDP(25) + (-g(1) * (-t69 * t104 + t70 * t72) - g(2) * (-t67 * t104 + t68 * t72) - (t74 * t101 + t72 * t82) * t106) * MDP(31) + (-g(1) * (t69 * t105 + t70 * t74) - g(2) * (t67 * t105 + t68 * t74) - (-t72 * t101 + t74 * t82) * t106) * MDP(32) + (-t84 * MDP(10) + t81 * MDP(11) - t75 * MDP(17) + t73 * MDP(18) - MDP(3)) * (-g(1) * t69 - g(2) * t67 + g(3) * t97); (-g(1) * (-t70 * t81 + t78 * t98) - g(2) * (-t68 * t81 - t84 * t90) - g(3) * (-t81 * t99 + t93 * t84)) * MDP(10) + (-g(1) * (-t81 * t100 - t70 * t84) - g(2) * (-t68 * t84 + t81 * t90) - g(3) * (-t93 * t81 - t82 * t98)) * MDP(11) + t88; t88; (-g(1) * (-t64 * t80 + t69 * t83) - g(2) * (-t62 * t80 + t67 * t83) - g(3) * (-t66 * t80 - t79 * t95)) * MDP(24) + (-g(1) * (-t64 * t83 - t69 * t80) - g(2) * (-t62 * t83 - t67 * t80) - g(3) * (-t66 * t83 + t79 * t96)) * MDP(25) + t94; t94;];
taug  = t1;
