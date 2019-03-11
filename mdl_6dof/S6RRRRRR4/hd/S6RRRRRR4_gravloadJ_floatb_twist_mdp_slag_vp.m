% Calculate Gravitation load on the joints for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:51:26
% EndTime: 2019-03-10 03:51:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (448->74), mult. (392->111), div. (0->0), fcn. (406->12), ass. (0->48)
t85 = sin(qJ(2));
t114 = t85 * MDP(10);
t83 = qJ(3) + qJ(4);
t82 = qJ(5) + t83;
t79 = qJ(6) + t82;
t75 = sin(t79);
t76 = cos(t79);
t77 = sin(t82);
t78 = cos(t82);
t80 = sin(t83);
t81 = cos(t83);
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t113 = t87 * MDP(16) - t84 * MDP(17) + t81 * MDP(23) - t80 * MDP(24) + t78 * MDP(30) - t77 * MDP(31) + t76 * MDP(37) - t75 * MDP(38) + MDP(9);
t112 = g(3) * t85;
t86 = sin(qJ(1));
t88 = cos(qJ(2));
t111 = t86 * t88;
t89 = cos(qJ(1));
t110 = t89 * t75;
t109 = t89 * t76;
t108 = t89 * t77;
t107 = t89 * t78;
t106 = t89 * t80;
t105 = t89 * t81;
t104 = t89 * t84;
t103 = t89 * t87;
t59 = t75 * t111 + t109;
t60 = -t76 * t111 + t110;
t61 = -t88 * t110 + t86 * t76;
t62 = t88 * t109 + t86 * t75;
t102 = (-g(1) * t61 + g(2) * t59 + t75 * t112) * MDP(37) + (g(1) * t62 - g(2) * t60 + t76 * t112) * MDP(38);
t63 = t77 * t111 + t107;
t64 = -t78 * t111 + t108;
t65 = -t88 * t108 + t86 * t78;
t66 = t88 * t107 + t86 * t77;
t93 = (-g(1) * t65 + g(2) * t63 + t77 * t112) * MDP(30) + (g(1) * t66 - g(2) * t64 + t78 * t112) * MDP(31) + t102;
t92 = g(1) * t89 + g(2) * t86;
t67 = t80 * t111 + t105;
t68 = -t81 * t111 + t106;
t69 = -t88 * t106 + t86 * t81;
t70 = t88 * t105 + t86 * t80;
t90 = (-g(1) * t69 + g(2) * t67 + t80 * t112) * MDP(23) + (g(1) * t70 - g(2) * t68 + t81 * t112) * MDP(24) + t93;
t74 = t88 * t103 + t86 * t84;
t73 = -t88 * t104 + t86 * t87;
t72 = -t87 * t111 + t104;
t71 = t84 * t111 + t103;
t1 = [t92 * MDP(3) + (-g(1) * t72 - g(2) * t74) * MDP(16) + (-g(1) * t71 - g(2) * t73) * MDP(17) + (-g(1) * t68 - g(2) * t70) * MDP(23) + (-g(1) * t67 - g(2) * t69) * MDP(24) + (-g(1) * t64 - g(2) * t66) * MDP(30) + (-g(1) * t63 - g(2) * t65) * MDP(31) + (-g(1) * t60 - g(2) * t62) * MDP(37) + (-g(1) * t59 - g(2) * t61) * MDP(38) + (t88 * MDP(9) + MDP(2) - t114) * (g(1) * t86 - g(2) * t89); (-t113 * t88 + t114) * g(3) + (MDP(10) * t88 + t113 * t85) * t92; (-g(1) * t73 + g(2) * t71 + t84 * t112) * MDP(16) + (g(1) * t74 - g(2) * t72 + t87 * t112) * MDP(17) + t90; t90; t93; t102;];
taug  = t1;
