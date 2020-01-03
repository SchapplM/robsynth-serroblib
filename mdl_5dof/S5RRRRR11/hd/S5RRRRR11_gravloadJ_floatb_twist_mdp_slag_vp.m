% Calculate Gravitation load on the joints for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:43
% EndTime: 2019-12-31 22:43:45
% DurationCPUTime: 0.48s
% Computational Cost: add. (295->93), mult. (630->165), div. (0->0), fcn. (770->12), ass. (0->44)
t81 = sin(qJ(2));
t82 = sin(qJ(1));
t85 = cos(qJ(2));
t106 = cos(qJ(1));
t92 = cos(pkin(5));
t89 = t92 * t106;
t70 = t81 * t89 + t82 * t85;
t80 = sin(qJ(3));
t84 = cos(qJ(3));
t78 = sin(pkin(5));
t91 = t78 * t106;
t63 = t70 * t84 - t80 * t91;
t69 = t82 * t81 - t85 * t89;
t79 = sin(qJ(4));
t83 = cos(qJ(4));
t111 = t63 * t79 - t69 * t83;
t110 = t63 * t83 + t69 * t79;
t77 = qJ(4) + qJ(5);
t75 = sin(t77);
t76 = cos(t77);
t109 = t63 * t75 - t69 * t76;
t108 = t63 * t76 + t69 * t75;
t107 = g(3) * t78;
t101 = t75 * t84;
t100 = t76 * t84;
t99 = t78 * t81;
t98 = t78 * t84;
t97 = t78 * t85;
t96 = t79 * t84;
t95 = t83 * t84;
t94 = t84 * t85;
t90 = t82 * t92;
t72 = t106 * t85 - t81 * t90;
t66 = t82 * t78 * t80 + t72 * t84;
t71 = t106 * t81 + t85 * t90;
t58 = -t66 * t75 + t71 * t76;
t59 = t66 * t76 + t71 * t75;
t68 = t80 * t92 + t81 * t98;
t93 = (-g(1) * t58 + g(2) * t109 - g(3) * (-t68 * t75 - t76 * t97)) * MDP(30) + (g(1) * t59 + g(2) * t108 - g(3) * (-t68 * t76 + t75 * t97)) * MDP(31);
t87 = t70 * t80 + t84 * t91;
t65 = -t72 * t80 + t82 * t98;
t61 = t66 * t83 + t71 * t79;
t60 = -t66 * t79 + t71 * t83;
t1 = [(g(1) * t82 - g(2) * t106) * MDP(2) + (g(1) * t106 + g(2) * t82) * MDP(3) + (g(1) * t70 - g(2) * t72) * MDP(9) + (-g(1) * t69 + g(2) * t71) * MDP(10) + (g(1) * t63 - g(2) * t66) * MDP(16) + (-g(1) * t87 - g(2) * t65) * MDP(17) + (g(1) * t110 - g(2) * t61) * MDP(23) + (-g(1) * t111 - g(2) * t60) * MDP(24) + (g(1) * t108 - g(2) * t59) * MDP(30) + (-g(1) * t109 - g(2) * t58) * MDP(31); (g(1) * t72 + g(2) * t70 + g(3) * t99) * MDP(10) + (-g(1) * (-t71 * t95 + t72 * t79) - g(2) * (-t69 * t95 + t70 * t79) - (t79 * t81 + t83 * t94) * t107) * MDP(23) + (-g(1) * (t71 * t96 + t72 * t83) - g(2) * (t69 * t96 + t70 * t83) - (-t79 * t94 + t81 * t83) * t107) * MDP(24) + (-g(1) * (-t100 * t71 + t72 * t75) - g(2) * (-t69 * t100 + t70 * t75) - (t75 * t81 + t76 * t94) * t107) * MDP(30) + (-g(1) * (t101 * t71 + t72 * t76) - g(2) * (t69 * t101 + t70 * t76) - (-t75 * t94 + t76 * t81) * t107) * MDP(31) + (-t84 * MDP(16) + t80 * MDP(17) - MDP(9)) * (-g(1) * t71 - g(2) * t69 + g(3) * t97); (g(1) * t66 + g(2) * t63 + g(3) * t68) * MDP(17) + (-MDP(23) * t83 + MDP(24) * t79 - MDP(30) * t76 + MDP(31) * t75 - MDP(16)) * (g(1) * t65 - g(2) * t87 + g(3) * (-t80 * t99 + t84 * t92)); (-g(1) * t60 + g(2) * t111 - g(3) * (-t68 * t79 - t83 * t97)) * MDP(23) + (g(1) * t61 + g(2) * t110 - g(3) * (-t68 * t83 + t79 * t97)) * MDP(24) + t93; t93;];
taug = t1;
