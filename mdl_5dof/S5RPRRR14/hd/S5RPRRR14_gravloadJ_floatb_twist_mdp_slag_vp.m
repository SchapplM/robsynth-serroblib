% Calculate Gravitation load on the joints for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:44
% EndTime: 2019-12-31 19:19:46
% DurationCPUTime: 0.52s
% Computational Cost: add. (376->87), mult. (1038->160), div. (0->0), fcn. (1327->14), ass. (0->53)
t107 = sin(qJ(5));
t111 = cos(qJ(5));
t108 = sin(qJ(4));
t112 = cos(qJ(4));
t109 = sin(qJ(3));
t135 = cos(pkin(6));
t106 = cos(pkin(5));
t113 = cos(qJ(1));
t134 = cos(pkin(11));
t122 = t113 * t134;
t103 = sin(pkin(11));
t110 = sin(qJ(1));
t129 = t110 * t103;
t95 = -t106 * t122 + t129;
t125 = t95 * t135;
t105 = sin(pkin(5));
t131 = t105 * t113;
t104 = sin(pkin(6));
t133 = t104 * t109;
t136 = cos(qJ(3));
t123 = t110 * t134;
t127 = t113 * t103;
t96 = t106 * t127 + t123;
t82 = -t109 * t125 - t131 * t133 + t96 * t136;
t124 = t105 * t135;
t90 = t95 * t104 - t113 * t124;
t75 = t90 * t108 + t82 * t112;
t126 = t105 * t136;
t121 = t104 * t126;
t81 = t96 * t109 + t113 * t121 + t136 * t125;
t141 = t75 * t107 - t81 * t111;
t140 = t81 * t107 + t75 * t111;
t137 = t82 * t108 - t90 * t112;
t132 = t105 * t110;
t130 = t107 * t112;
t128 = t111 * t112;
t119 = g(1) * t110 - g(2) * t113;
t118 = t135 * t134;
t115 = t106 * t123 + t127;
t114 = t115 * t135;
t97 = -t106 * t129 + t122;
t94 = -t105 * t134 * t104 + t106 * t135;
t91 = t115 * t104 + t110 * t124;
t88 = t106 * t133 + (t136 * t103 + t109 * t118) * t105;
t87 = t105 * t103 * t109 - t106 * t104 * t136 - t118 * t126;
t86 = t97 * t136 + (t104 * t132 - t114) * t109;
t85 = t97 * t109 - t110 * t121 + t136 * t114;
t80 = t94 * t108 + t88 * t112;
t78 = t91 * t108 + t86 * t112;
t77 = -t86 * t108 + t91 * t112;
t73 = t85 * t107 + t78 * t111;
t72 = -t78 * t107 + t85 * t111;
t1 = [t119 * MDP(2) + (g(1) * t96 - g(2) * t97) * MDP(4) + (-g(1) * t95 + g(2) * t115) * MDP(5) + (-g(1) * (-t110 * pkin(1) + qJ(2) * t131) - g(2) * (t113 * pkin(1) + qJ(2) * t132)) * MDP(7) + (g(1) * t82 - g(2) * t86) * MDP(13) + (-g(1) * t81 + g(2) * t85) * MDP(14) + (g(1) * t75 - g(2) * t78) * MDP(20) + (-g(1) * t137 - g(2) * t77) * MDP(21) + (g(1) * t140 - g(2) * t73) * MDP(27) + (-g(1) * t141 - g(2) * t72) * MDP(28) + (-t105 * MDP(6) + MDP(3)) * (g(1) * t113 + g(2) * t110); (-g(3) * t106 - t119 * t105) * MDP(7); (g(1) * t86 + g(2) * t82 + g(3) * t88) * MDP(14) + (-g(1) * (t86 * t107 - t85 * t128) - g(2) * (t82 * t107 - t81 * t128) - g(3) * (t88 * t107 - t87 * t128)) * MDP(27) + (-g(1) * (t86 * t111 + t85 * t130) - g(2) * (t82 * t111 + t81 * t130) - g(3) * (t88 * t111 + t87 * t130)) * MDP(28) + (t112 * MDP(20) - MDP(21) * t108 + MDP(13)) * (g(1) * t85 + g(2) * t81 + g(3) * t87); (g(1) * t78 + g(2) * t75 + g(3) * t80) * MDP(21) + (-MDP(27) * t111 + MDP(28) * t107 - MDP(20)) * (g(1) * t77 - g(2) * t137 + g(3) * (-t88 * t108 + t94 * t112)); (-g(1) * t72 + g(2) * t141 - g(3) * (-t80 * t107 + t87 * t111)) * MDP(27) + (g(1) * t73 + g(2) * t140 - g(3) * (-t87 * t107 - t80 * t111)) * MDP(28);];
taug = t1;
