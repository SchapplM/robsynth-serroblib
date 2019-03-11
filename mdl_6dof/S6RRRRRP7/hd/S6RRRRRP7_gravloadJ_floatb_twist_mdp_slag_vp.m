% Calculate Gravitation load on the joints for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:11
% EndTime: 2019-03-10 01:43:13
% DurationCPUTime: 0.72s
% Computational Cost: add. (538->119), mult. (905->189), div. (0->0), fcn. (1057->12), ass. (0->59)
t176 = MDP(24) - MDP(32);
t130 = sin(qJ(2));
t131 = sin(qJ(1));
t134 = cos(qJ(2));
t135 = cos(qJ(1));
t163 = cos(pkin(6));
t148 = t135 * t163;
t112 = t130 * t148 + t131 * t134;
t125 = qJ(3) + qJ(4);
t123 = sin(t125);
t124 = cos(t125);
t126 = sin(pkin(6));
t154 = t126 * t135;
t100 = t112 * t124 - t123 * t154;
t111 = t130 * t131 - t134 * t148;
t128 = sin(qJ(5));
t132 = cos(qJ(5));
t175 = t100 * t128 - t111 * t132;
t162 = t111 * t128;
t174 = t100 * t132 + t162;
t149 = t131 * t163;
t114 = -t130 * t149 + t134 * t135;
t129 = sin(qJ(3));
t133 = cos(qJ(3));
t155 = t126 * t133;
t105 = -t114 * t129 + t131 * t155;
t144 = t112 * t129 + t133 * t154;
t157 = t126 * t130;
t173 = g(2) * t144 - g(3) * (-t129 * t157 + t163 * t133) - g(1) * t105;
t171 = g(1) * t135 + g(2) * t131;
t167 = g(2) * t111;
t166 = g(2) * t112;
t164 = g(3) * t126;
t113 = t135 * t130 + t134 * t149;
t160 = t113 * t128;
t159 = t124 * t128;
t158 = t124 * t132;
t156 = t126 * t131;
t153 = t128 * t134;
t152 = t132 * t134;
t136 = -pkin(10) - pkin(9);
t150 = pkin(5) * t128 - t136;
t147 = t112 * t133 - t129 * t154;
t104 = t114 * t124 + t123 * t156;
t110 = t163 * t123 + t124 * t157;
t103 = t114 * t123 - t124 * t156;
t109 = t123 * t157 - t163 * t124;
t99 = t112 * t123 + t124 * t154;
t142 = g(1) * t103 + g(2) * t99 + g(3) * t109;
t146 = t176 * (g(1) * t104 + g(2) * t100 + g(3) * t110) + (MDP(30) * t132 - MDP(31) * t128 + MDP(23)) * t142;
t93 = -t104 * t128 + t113 * t132;
t121 = pkin(5) * t132 + pkin(4);
t122 = pkin(3) * t133 + pkin(2);
t127 = -qJ(6) - pkin(11);
t143 = t121 * t124 - t123 * t127 + t122;
t139 = -g(1) * (-t103 * t121 - t104 * t127) - g(2) * (-t100 * t127 - t99 * t121) - g(3) * (-t109 * t121 - t110 * t127);
t106 = t114 * t133 + t129 * t156;
t94 = t104 * t132 + t160;
t1 = [(g(1) * t131 - g(2) * t135) * MDP(2) + t171 * MDP(3) + (g(1) * t112 - g(2) * t114) * MDP(9) + (-g(1) * t111 + g(2) * t113) * MDP(10) + (g(1) * t147 - g(2) * t106) * MDP(16) + (-g(1) * t144 - g(2) * t105) * MDP(17) + (g(1) * t100 - g(2) * t104) * MDP(23) + (g(1) * t174 - g(2) * t94) * MDP(30) + (-g(1) * t175 - g(2) * t93) * MDP(31) + (-g(1) * (-t131 * pkin(1) - pkin(5) * t162 - t100 * t121 + t111 * t136 - t112 * t122 + t127 * t99) - g(2) * (t135 * pkin(1) + pkin(5) * t160 - t103 * t127 + t104 * t121 - t113 * t136 + t114 * t122) - t171 * t126 * (pkin(3) * t129 + pkin(8))) * MDP(33) + t176 * (-g(1) * t99 + g(2) * t103); (g(1) * t114 + g(3) * t157 + t166) * MDP(10) + (-g(1) * (-t113 * t158 + t114 * t128) - g(2) * (-t111 * t158 + t112 * t128) - (t124 * t152 + t128 * t130) * t164) * MDP(30) + (-g(1) * (t113 * t159 + t114 * t132) - g(2) * (t111 * t159 + t112 * t132) - (-t124 * t153 + t130 * t132) * t164) * MDP(31) + (-g(1) * (-t143 * t113 + t114 * t150) - t150 * t166 + t143 * t167 - (t130 * t150 + t143 * t134) * t164) * MDP(33) + (-t133 * MDP(16) + t129 * MDP(17) - t124 * MDP(23) + t176 * t123 - MDP(9)) * (-g(1) * t113 + t134 * t164 - t167); t173 * MDP(16) + (g(1) * t106 + g(2) * t147 - g(3) * (-t163 * t129 - t130 * t155)) * MDP(17) + (pkin(3) * t173 + t139) * MDP(33) + t146; MDP(33) * t139 + t146; (g(1) * t94 + g(2) * t174 - g(3) * (-t110 * t132 + t126 * t153)) * MDP(31) + (MDP(33) * pkin(5) + MDP(30)) * (g(2) * t175 - g(3) * (-t110 * t128 - t126 * t152) - g(1) * t93); -t142 * MDP(33);];
taug  = t1;
