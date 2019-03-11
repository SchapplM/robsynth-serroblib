% Calculate Gravitation load on the joints for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:48
% EndTime: 2019-03-09 20:46:50
% DurationCPUTime: 0.52s
% Computational Cost: add. (451->101), mult. (495->135), div. (0->0), fcn. (454->10), ass. (0->56)
t116 = qJ(4) + pkin(10);
t111 = sin(t116);
t112 = cos(t116);
t160 = pkin(5) * t112 + qJ(6) * t111;
t123 = cos(qJ(2));
t115 = t123 * pkin(2);
t117 = qJ(2) + qJ(3);
t113 = sin(t117);
t118 = -qJ(5) - pkin(9);
t148 = t113 * t118;
t159 = t115 + pkin(1) - t148;
t158 = MDP(17) - MDP(25) - MDP(28);
t114 = cos(t117);
t157 = t160 * t114;
t121 = sin(qJ(1));
t124 = cos(qJ(1));
t137 = g(1) * t124 + g(2) * t121;
t86 = -g(3) * t114 + t137 * t113;
t120 = sin(qJ(2));
t156 = pkin(2) * t120;
t152 = g(3) * t113;
t149 = t111 * t121;
t122 = cos(qJ(4));
t109 = pkin(4) * t122 + pkin(3);
t100 = t114 * t109;
t147 = t114 * t118;
t146 = t114 * t124;
t119 = sin(qJ(4));
t145 = t119 * t121;
t144 = t119 * t124;
t143 = t121 * t112;
t142 = t121 * t122;
t141 = t122 * t124;
t140 = t124 * t111;
t139 = t114 * t144;
t138 = t100 - t148;
t135 = t115 + t138;
t134 = t109 * t113 + t147;
t133 = t109 + t160;
t132 = t137 * t114;
t93 = t114 * t145 + t141;
t87 = t132 + t152;
t131 = t158 * t87 + (MDP(23) * t122 - MDP(24) * t119 + MDP(27) * t112 + MDP(29) * t111 + MDP(16)) * t86;
t88 = t112 * t124 + t114 * t149;
t90 = t114 * t140 - t143;
t129 = g(1) * t90 + g(2) * t88 + t111 * t152;
t125 = -pkin(8) - pkin(7);
t128 = pkin(4) * t145 + t109 * t146 - t121 * t125 + t159 * t124;
t127 = pkin(4) * t144 - t124 * t125 + (-t100 - t159) * t121;
t106 = pkin(4) * t142;
t96 = t114 * t141 + t145;
t95 = -t139 + t142;
t94 = -t114 * t142 + t144;
t91 = t112 * t146 + t149;
t89 = t114 * t143 - t140;
t1 = [t137 * MDP(3) + (-g(1) * t94 - g(2) * t96) * MDP(23) + (-g(1) * t93 - g(2) * t95) * MDP(24) + (-g(1) * t127 - g(2) * t128) * MDP(26) + (g(1) * t89 - g(2) * t91) * MDP(27) + (g(1) * t88 - g(2) * t90) * MDP(29) + (-g(1) * (-t89 * pkin(5) - t88 * qJ(6) + t127) - g(2) * (t91 * pkin(5) + t90 * qJ(6) + t128)) * MDP(30) + (-t120 * MDP(10) + t114 * MDP(16) + t123 * MDP(9) - t158 * t113 + MDP(2)) * (g(1) * t121 - g(2) * t124); (-g(3) * t123 + t137 * t120) * MDP(9) + (g(3) * t120 + t137 * t123) * MDP(10) + (-g(3) * t135 + t137 * (t134 + t156)) * MDP(26) + (-g(3) * (t135 + t157) + t137 * (t133 * t113 + t147 + t156)) * MDP(30) + t131; (-g(3) * t138 + t137 * t134) * MDP(26) + (-g(3) * (t100 + t157) + t118 * t132 + (g(3) * t118 + t137 * t133) * t113) * MDP(30) + t131; (-g(1) * t95 + g(2) * t93 + t119 * t152) * MDP(23) + (g(1) * t96 - g(2) * t94 + t122 * t152) * MDP(24) + (-g(1) * t106 + (g(2) * t141 + t119 * t87) * pkin(4)) * MDP(26) + t129 * MDP(27) + (-g(1) * t91 - g(2) * t89 - t112 * t152) * MDP(29) + (-g(1) * (-pkin(4) * t139 - pkin(5) * t90 + qJ(6) * t91 + t106) - g(2) * (-pkin(4) * t93 - pkin(5) * t88 + qJ(6) * t89) - (-pkin(4) * t119 - pkin(5) * t111 + qJ(6) * t112) * t152) * MDP(30); (-MDP(26) - MDP(30)) * t86; -t129 * MDP(30);];
taug  = t1;
