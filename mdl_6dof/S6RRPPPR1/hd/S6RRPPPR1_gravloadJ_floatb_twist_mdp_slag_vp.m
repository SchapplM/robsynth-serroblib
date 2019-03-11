% Calculate Gravitation load on the joints for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:19
% EndTime: 2019-03-09 08:08:21
% DurationCPUTime: 0.60s
% Computational Cost: add. (281->87), mult. (405->119), div. (0->0), fcn. (399->10), ass. (0->51)
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t83 = g(1) * t104 + g(2) * t101;
t95 = qJ(2) + pkin(9);
t91 = sin(t95);
t131 = t83 * t91;
t92 = cos(t95);
t107 = -g(3) * t92 + t131;
t135 = MDP(13) + MDP(17);
t134 = -MDP(14) + MDP(19);
t133 = MDP(15) + MDP(18);
t87 = t92 * pkin(3);
t132 = -t91 * qJ(4) - t87;
t129 = g(3) * t91;
t100 = sin(qJ(2));
t127 = pkin(2) * t100;
t96 = sin(pkin(10));
t124 = t101 * t96;
t97 = cos(pkin(10));
t123 = t101 * t97;
t122 = t104 * t96;
t121 = t104 * t97;
t98 = -qJ(3) - pkin(7);
t120 = t104 * t98;
t119 = qJ(4) * t104;
t118 = MDP(16) + MDP(20);
t103 = cos(qJ(2));
t93 = t103 * pkin(2);
t90 = t93 + pkin(1);
t85 = t104 * t90;
t117 = t104 * t87 + t91 * t119 + t85;
t116 = t93 - t132;
t114 = -pkin(3) * t91 - t127;
t113 = pkin(4) * t97 + qJ(5) * t96;
t82 = g(1) * t101 - g(2) * t104;
t102 = cos(qJ(6));
t75 = t92 * t124 + t121;
t76 = t92 * t123 - t122;
t99 = sin(qJ(6));
t112 = t102 * t76 + t75 * t99;
t111 = t102 * t75 - t76 * t99;
t110 = t102 * t97 + t96 * t99;
t109 = t102 * t96 - t97 * t99;
t105 = (-g(1) * (-t90 + t132) + g(2) * t98) * t101;
t81 = t92 * t119;
t79 = t101 * t92 * qJ(4);
t78 = t92 * t121 + t124;
t77 = t92 * t122 - t123;
t68 = t102 * t78 + t77 * t99;
t67 = t102 * t77 - t78 * t99;
t1 = [(-g(1) * (-t101 * t90 - t120) - g(2) * (-t101 * t98 + t85)) * MDP(12) + (g(1) * t120 - g(2) * t117 + t105) * MDP(16) + (-g(1) * (-t76 * pkin(4) - t75 * qJ(5) - t120) - g(2) * (pkin(4) * t78 + qJ(5) * t77 + t117) + t105) * MDP(20) + (g(1) * t112 - g(2) * t68) * MDP(26) + (g(1) * t111 - g(2) * t67) * MDP(27) + (MDP(3) - MDP(11)) * t83 + t135 * (g(1) * t76 - g(2) * t78) + t134 * (g(1) * t75 - g(2) * t77) + (-t100 * MDP(10) + t103 * MDP(9) + t133 * t91 + MDP(2)) * t82; (g(3) * t100 + t83 * t103) * MDP(10) + (-g(1) * (t114 * t104 + t81) - g(2) * (t114 * t101 + t79) - g(3) * t116) * MDP(16) + (-g(1) * (-t104 * t127 + t81) - g(2) * (-t101 * t127 + t79) - g(3) * (t113 * t92 + t116) + (pkin(3) + t113) * t131) * MDP(20) + t133 * (-t83 * t92 - t129) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t103 + t83 * t100) + (t110 * MDP(26) + t109 * MDP(27) + t134 * t96 + t135 * t97) * t107; (-MDP(12) - t118) * t82; -t118 * t107; (-g(1) * t77 - g(2) * t75 - t96 * t129) * MDP(20); (-g(1) * t67 - g(2) * t111) * MDP(26) + (g(1) * t68 + g(2) * t112) * MDP(27) + (-t109 * MDP(26) + t110 * MDP(27)) * t129;];
taug  = t1;
