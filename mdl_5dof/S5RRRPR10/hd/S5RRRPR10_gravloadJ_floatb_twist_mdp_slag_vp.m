% Calculate Gravitation load on the joints for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:52
% EndTime: 2019-12-31 21:29:54
% DurationCPUTime: 0.56s
% Computational Cost: add. (231->89), mult. (470->156), div. (0->0), fcn. (542->12), ass. (0->45)
t125 = MDP(10) - MDP(18);
t85 = sin(pkin(5));
t94 = cos(qJ(1));
t109 = t85 * t94;
t106 = cos(pkin(5));
t102 = t94 * t106;
t89 = sin(qJ(2));
t90 = sin(qJ(1));
t93 = cos(qJ(2));
t75 = t89 * t102 + t90 * t93;
t84 = qJ(3) + pkin(10);
t82 = sin(t84);
t83 = cos(t84);
t67 = t82 * t109 - t75 * t83;
t74 = -t93 * t102 + t90 * t89;
t87 = sin(qJ(5));
t91 = cos(qJ(5));
t124 = t67 * t87 + t74 * t91;
t123 = t67 * t91 - t74 * t87;
t122 = g(1) * t94 + g(2) * t90;
t117 = g(3) * t85;
t114 = t83 * t87;
t113 = t83 * t91;
t112 = t85 * t89;
t111 = t85 * t90;
t92 = cos(qJ(3));
t110 = t85 * t92;
t108 = t87 * t93;
t107 = t91 * t93;
t88 = sin(qJ(3));
t104 = -t88 * t109 + t75 * t92;
t103 = t90 * t106;
t100 = t92 * t109 + t75 * t88;
t77 = -t89 * t103 + t94 * t93;
t70 = t90 * t110 - t77 * t88;
t76 = t93 * t103 + t94 * t89;
t97 = -g(1) * t76 - g(2) * t74 + t93 * t117;
t86 = -qJ(4) - pkin(8);
t81 = t92 * pkin(3) + pkin(2);
t73 = t106 * t82 + t83 * t112;
t71 = t88 * t111 + t77 * t92;
t69 = t82 * t111 + t77 * t83;
t64 = t69 * t91 + t76 * t87;
t63 = -t69 * t87 + t76 * t91;
t1 = [(g(1) * t90 - g(2) * t94) * MDP(2) + t122 * MDP(3) + (g(1) * t75 - g(2) * t77) * MDP(9) + (g(1) * t104 - g(2) * t71) * MDP(16) + (-g(1) * t100 - g(2) * t70) * MDP(17) + (-g(1) * (-t90 * pkin(1) + t74 * t86 - t75 * t81) - g(2) * (t94 * pkin(1) - t76 * t86 + t77 * t81) - t122 * t85 * (pkin(3) * t88 + pkin(7))) * MDP(19) + (-g(1) * t123 - g(2) * t64) * MDP(25) + (g(1) * t124 - g(2) * t63) * MDP(26) - t125 * (g(1) * t74 - g(2) * t76); (-g(1) * (-t76 * t81 - t77 * t86) - g(2) * (-t74 * t81 - t75 * t86) - (t81 * t93 - t86 * t89) * t117) * MDP(19) + (-g(1) * (-t76 * t113 + t77 * t87) - g(2) * (-t74 * t113 + t75 * t87) - (t83 * t107 + t87 * t89) * t117) * MDP(25) + (-g(1) * (t76 * t114 + t77 * t91) - g(2) * (t74 * t114 + t75 * t91) - (-t83 * t108 + t89 * t91) * t117) * MDP(26) + t125 * (g(1) * t77 + g(2) * t75 + g(3) * t112) + (-t92 * MDP(16) + t88 * MDP(17) - MDP(9)) * t97; (g(1) * t71 + g(2) * t104 - g(3) * (-t106 * t88 - t89 * t110)) * MDP(17) + (-MDP(25) * t91 + MDP(26) * t87) * (g(1) * (t83 * t111 - t77 * t82) + g(2) * (-t83 * t109 - t75 * t82) + g(3) * (t106 * t83 - t82 * t112)) + (pkin(3) * MDP(19) + MDP(16)) * (-g(3) * (t106 * t92 - t88 * t112) + g(2) * t100 - g(1) * t70); t97 * MDP(19); (-g(1) * t63 - g(2) * t124 - g(3) * (-t85 * t107 - t73 * t87)) * MDP(25) + (g(1) * t64 - g(2) * t123 - g(3) * (t85 * t108 - t73 * t91)) * MDP(26);];
taug = t1;
