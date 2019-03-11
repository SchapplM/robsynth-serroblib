% Calculate Gravitation load on the joints for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:37
% EndTime: 2019-03-09 12:48:38
% DurationCPUTime: 0.40s
% Computational Cost: add. (217->86), mult. (349->121), div. (0->0), fcn. (317->8), ass. (0->45)
t105 = sin(qJ(2));
t108 = cos(qJ(2));
t123 = t108 * pkin(2) + t105 * qJ(3);
t106 = sin(qJ(1));
t109 = cos(qJ(1));
t87 = g(1) * t109 + g(2) * t106;
t135 = MDP(10) - MDP(13);
t126 = g(3) * t108;
t120 = t105 * t109;
t103 = qJ(4) + qJ(5);
t93 = sin(t103);
t94 = cos(t103);
t72 = -t106 * t93 + t94 * t120;
t121 = t105 * t106;
t74 = t109 * t93 + t94 * t121;
t134 = -g(1) * t72 - g(2) * t74 + t94 * t126;
t133 = MDP(9) - MDP(12) + MDP(29);
t77 = g(3) * t105 + t87 * t108;
t130 = g(1) * t106;
t73 = t106 * t94 + t93 * t120;
t75 = t109 * t94 - t93 * t121;
t124 = t134 * MDP(27) + (g(1) * t73 - g(2) * t75 - t93 * t126) * MDP(28);
t107 = cos(qJ(4));
t86 = t107 * pkin(4) + pkin(5) * t94;
t102 = -qJ(6) - pkin(9) - pkin(8);
t122 = t102 * t108;
t104 = sin(qJ(4));
t119 = t106 * t104;
t118 = t106 * t107;
t117 = t106 * t108;
t116 = t107 * t109;
t115 = t108 * t109;
t113 = t109 * pkin(1) + pkin(2) * t115 + t106 * pkin(7) + qJ(3) * t120;
t111 = -pkin(1) - t123;
t99 = t109 * pkin(7);
t90 = qJ(3) * t115;
t88 = qJ(3) * t117;
t85 = pkin(4) * t104 + pkin(5) * t93;
t84 = pkin(3) + t86;
t81 = -t105 * t119 + t116;
t80 = t104 * t109 + t105 * t118;
t79 = t104 * t120 + t118;
t78 = t105 * t116 - t119;
t76 = t87 * t105 - t126;
t1 = [(-g(1) * t99 - g(2) * t113 - t111 * t130) * MDP(14) + (-g(1) * t81 - g(2) * t79) * MDP(20) + (g(1) * t80 - g(2) * t78) * MDP(21) + (-g(1) * t75 - g(2) * t73) * MDP(27) + (g(1) * t74 - g(2) * t72) * MDP(28) + (-g(1) * (t109 * t84 + t99) - g(2) * (-t102 * t115 + t85 * t120 + t113) + (-g(1) * (-t105 * t85 + t111 + t122) - g(2) * t84) * t106) * MDP(30) + (MDP(3) - MDP(11)) * t87 + (-t135 * t105 + t133 * t108 + MDP(2)) * (-g(2) * t109 + t130); (-g(1) * (-pkin(2) * t120 + t90) - g(2) * (-pkin(2) * t121 + t88) - g(3) * t123) * MDP(14) + (-g(1) * (t85 * t115 + t90) - g(2) * (t85 * t117 + t88) - g(3) * (-t122 + t123) + (-g(3) * t85 + t87 * (pkin(2) - t102)) * t105) * MDP(30) + t133 * t76 + (-t104 * MDP(20) - t107 * MDP(21) - t93 * MDP(27) - t94 * MDP(28) + t135) * t77; (-MDP(14) - MDP(30)) * t76; (-g(1) * t78 - g(2) * t80 + t107 * t126) * MDP(20) + (g(1) * t79 - g(2) * t81 - t104 * t126) * MDP(21) + (-g(1) * (-t106 * t85 + t86 * t120) - g(2) * (t109 * t85 + t86 * t121) + t86 * t126) * MDP(30) + t124; t134 * MDP(30) * pkin(5) + t124; -t77 * MDP(30);];
taug  = t1;
