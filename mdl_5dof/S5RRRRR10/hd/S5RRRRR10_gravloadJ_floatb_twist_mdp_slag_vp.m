% Calculate Gravitation load on the joints for
% S5RRRRR10
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
%   see S5RRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:59
% EndTime: 2019-12-31 22:36:00
% DurationCPUTime: 0.37s
% Computational Cost: add. (295->81), mult. (532->143), div. (0->0), fcn. (634->12), ass. (0->43)
t81 = sin(pkin(5));
t89 = cos(qJ(1));
t101 = t81 * t89;
t84 = sin(qJ(2));
t85 = sin(qJ(1));
t88 = cos(qJ(2));
t98 = cos(pkin(5));
t95 = t89 * t98;
t72 = t84 * t95 + t85 * t88;
t80 = qJ(3) + qJ(4);
t78 = sin(t80);
t79 = cos(t80);
t63 = -t78 * t101 + t72 * t79;
t71 = t85 * t84 - t88 * t95;
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t111 = t63 * t82 - t71 * t86;
t110 = t63 * t86 + t71 * t82;
t109 = g(3) * t81;
t106 = t79 * t82;
t105 = t79 * t86;
t104 = t81 * t84;
t103 = t81 * t85;
t87 = cos(qJ(3));
t102 = t81 * t87;
t100 = t82 * t88;
t99 = t86 * t88;
t83 = sin(qJ(3));
t97 = -t83 * t101 + t72 * t87;
t96 = t85 * t98;
t74 = -t84 * t96 + t89 * t88;
t65 = t79 * t103 - t74 * t78;
t66 = t78 * t103 + t74 * t79;
t70 = t79 * t104 + t98 * t78;
t93 = t79 * t101 + t72 * t78;
t94 = (g(1) * t66 + g(2) * t63 + g(3) * t70) * MDP(24) + (-MDP(30) * t86 + MDP(31) * t82 - MDP(23)) * (g(1) * t65 - g(2) * t93 + g(3) * (-t78 * t104 + t98 * t79));
t92 = t87 * t101 + t72 * t83;
t73 = t89 * t84 + t88 * t96;
t68 = t83 * t103 + t74 * t87;
t67 = t85 * t102 - t74 * t83;
t61 = t66 * t86 + t73 * t82;
t60 = -t66 * t82 + t73 * t86;
t1 = [(g(1) * t85 - g(2) * t89) * MDP(2) + (g(1) * t89 + g(2) * t85) * MDP(3) + (g(1) * t72 - g(2) * t74) * MDP(9) + (-g(1) * t71 + g(2) * t73) * MDP(10) + (g(1) * t97 - g(2) * t68) * MDP(16) + (-g(1) * t92 - g(2) * t67) * MDP(17) + (g(1) * t63 - g(2) * t66) * MDP(23) + (-g(1) * t93 - g(2) * t65) * MDP(24) + (g(1) * t110 - g(2) * t61) * MDP(30) + (-g(1) * t111 - g(2) * t60) * MDP(31); (g(1) * t74 + g(2) * t72 + g(3) * t104) * MDP(10) + (-g(1) * (-t73 * t105 + t74 * t82) - g(2) * (-t71 * t105 + t72 * t82) - (t79 * t99 + t82 * t84) * t109) * MDP(30) + (-g(1) * (t73 * t106 + t74 * t86) - g(2) * (t71 * t106 + t72 * t86) - (-t79 * t100 + t84 * t86) * t109) * MDP(31) + (-MDP(16) * t87 + MDP(17) * t83 - t79 * MDP(23) + MDP(24) * t78 - MDP(9)) * (-g(1) * t73 - g(2) * t71 + t88 * t109); (-g(1) * t67 + g(2) * t92 - g(3) * (-t83 * t104 + t98 * t87)) * MDP(16) + (g(1) * t68 + g(2) * t97 - g(3) * (-t84 * t102 - t98 * t83)) * MDP(17) + t94; t94; (-g(1) * t60 + g(2) * t111 - g(3) * (-t70 * t82 - t81 * t99)) * MDP(30) + (g(1) * t61 + g(2) * t110 - g(3) * (t81 * t100 - t70 * t86)) * MDP(31);];
taug = t1;
