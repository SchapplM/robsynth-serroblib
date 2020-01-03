% Calculate Gravitation load on the joints for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:40
% EndTime: 2019-12-31 20:26:42
% DurationCPUTime: 0.46s
% Computational Cost: add. (251->80), mult. (640->145), div. (0->0), fcn. (795->12), ass. (0->48)
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t109 = sin(qJ(1));
t113 = cos(qJ(1));
t105 = cos(pkin(5));
t102 = sin(pkin(10));
t104 = cos(pkin(10));
t108 = sin(qJ(2));
t112 = cos(qJ(2));
t119 = t112 * t102 + t108 * t104;
t90 = t119 * t105;
t96 = t108 * t102 - t112 * t104;
t121 = -t109 * t96 + t113 * t90;
t103 = sin(pkin(5));
t133 = t103 * t113;
t75 = -t107 * t133 + t111 * t121;
t114 = t96 * t105;
t80 = -t109 * t119 - t113 * t114;
t141 = t75 * t106 + t80 * t110;
t140 = -t80 * t106 + t75 * t110;
t137 = g(3) * t103;
t134 = t103 * t109;
t132 = t106 * t111;
t130 = t109 * t108;
t129 = t109 * t112;
t128 = t110 * t111;
t126 = t113 * t108;
t125 = t113 * t112;
t122 = g(1) * t109 - g(2) * t113;
t120 = -t109 * t90 - t113 * t96;
t118 = t107 * t121 + t111 * t133;
t117 = t105 * t125 - t130;
t94 = -t105 * t129 - t126;
t101 = t112 * pkin(2) + pkin(1);
t95 = -t105 * t130 + t125;
t93 = -t105 * t126 - t129;
t91 = t105 * t108 * pkin(2) + (-pkin(7) - qJ(3)) * t103;
t89 = t119 * t103;
t88 = t96 * t103;
t86 = t105 * t107 + t89 * t111;
t83 = t109 * t114 - t113 * t119;
t78 = t107 * t134 + t111 * t120;
t77 = -t107 * t120 + t111 * t134;
t73 = -t83 * t106 + t78 * t110;
t72 = -t78 * t106 - t83 * t110;
t1 = [t122 * MDP(2) + (-g(1) * t93 - g(2) * t95) * MDP(9) + (g(1) * t117 - g(2) * t94) * MDP(10) + (-g(1) * (-t109 * t101 - t113 * t91) - g(2) * (t113 * t101 - t109 * t91)) * MDP(12) + (g(1) * t75 - g(2) * t78) * MDP(18) + (-g(1) * t118 - g(2) * t77) * MDP(19) + (g(1) * t140 - g(2) * t73) * MDP(25) + (-g(1) * t141 - g(2) * t72) * MDP(26) + (-t103 * MDP(11) + MDP(3)) * (g(1) * t113 + g(2) * t109); (g(1) * t95 - g(2) * t93 + t108 * t137) * MDP(10) + (-g(1) * (t106 * t120 + t83 * t128) - g(2) * (t106 * t121 + t80 * t128) - g(3) * (t89 * t106 - t88 * t128)) * MDP(25) + (-g(1) * (t110 * t120 - t83 * t132) - g(2) * (t110 * t121 - t80 * t132) - g(3) * (t89 * t110 + t88 * t132)) * MDP(26) + (-t111 * MDP(18) + MDP(19) * t107) * (g(1) * t83 + g(2) * t80 - g(3) * t88) + (pkin(2) * MDP(12) + MDP(9)) * (-g(1) * t94 - g(2) * t117 - t112 * t137); (-g(3) * t105 - t122 * t103) * MDP(12); (g(1) * t78 + g(2) * t75 + g(3) * t86) * MDP(19) + (-MDP(25) * t110 + MDP(26) * t106 - MDP(18)) * (g(1) * t77 - g(2) * t118 + g(3) * (t105 * t111 - t89 * t107)); (-g(1) * t72 + g(2) * t141 - g(3) * (-t86 * t106 + t88 * t110)) * MDP(25) + (g(1) * t73 + g(2) * t140 - g(3) * (-t88 * t106 - t86 * t110)) * MDP(26);];
taug = t1;
