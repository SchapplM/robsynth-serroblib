% Calculate Gravitation load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:52
% EndTime: 2021-01-16 01:51:57
% DurationCPUTime: 1.04s
% Computational Cost: add. (279->122), mult. (694->191), div. (0->0), fcn. (799->10), ass. (0->72)
t178 = MDP(14) - MDP(24);
t126 = sin(pkin(6));
t177 = g(3) * t126;
t128 = cos(pkin(6));
t134 = cos(qJ(4));
t131 = sin(qJ(4));
t135 = cos(qJ(2));
t152 = t131 * t135;
t116 = -t126 * t152 + t128 * t134;
t130 = sin(qJ(5));
t132 = sin(qJ(2));
t133 = cos(qJ(5));
t151 = t132 * t133;
t175 = g(3) * (-t116 * t130 + t126 * t151);
t162 = t126 * t132;
t174 = g(3) * (-t116 * t133 - t130 * t162);
t154 = t131 * t132;
t173 = (-t130 * t154 + t133 * t135) * t177;
t172 = (t130 * t135 + t131 * t151) * t177;
t149 = t134 * t135;
t115 = t126 * t149 + t128 * t131;
t171 = g(3) * t115;
t125 = sin(pkin(10));
t127 = cos(pkin(10));
t157 = t128 * t135;
t111 = t125 * t132 - t127 * t157;
t107 = t111 * t130;
t158 = t128 * t132;
t112 = t125 * t135 + t127 * t158;
t169 = t112 * t130;
t109 = t112 * t133;
t113 = t125 * t157 + t127 * t132;
t168 = t113 * t130;
t167 = t113 * t133;
t143 = t125 * t158 - t127 * t135;
t166 = t143 * t130;
t165 = t143 * t133;
t164 = t125 * t128;
t163 = t126 * t131;
t161 = t126 * t134;
t160 = t126 * t135;
t159 = t127 * t128;
t136 = pkin(2) + pkin(8);
t156 = t128 * t136;
t155 = t130 * t131;
t153 = t131 * t133;
t150 = t132 * t134;
t148 = MDP(25) + MDP(7);
t147 = qJ(3) * t159;
t142 = t128 * t149 - t163;
t146 = g(1) * (t142 * t125 + t127 * t150) - g(2) * (-t125 * t150 + t142 * t127);
t95 = -t113 * t134 + t125 * t163;
t97 = t111 * t134 + t127 * t163;
t145 = -g(1) * t95 + g(2) * t97;
t122 = pkin(5) * t133 + pkin(4);
t129 = -qJ(6) - pkin(9);
t144 = t131 * t122 + t134 * t129;
t141 = t128 * t152 + t161;
t139 = -t145 + t171;
t90 = -g(1) * t113 - g(2) * t111 + g(3) * t160;
t124 = t127 * qJ(3);
t123 = t125 * qJ(3);
t121 = qJ(3) * t162;
t120 = qJ(3) * t164;
t108 = t111 * t133;
t102 = t112 * t131;
t101 = t143 * t131;
t98 = -t111 * t131 + t127 * t161;
t96 = t113 * t131 + t125 * t161;
t93 = -t125 * t154 + t141 * t127;
t91 = t141 * t125 + t127 * t154;
t1 = [(-MDP(1) - t148) * g(3); (-g(1) * (-(pkin(2) * t127 + t120) * t132 + (-pkin(2) * t164 + t124) * t135) - g(2) * (-(pkin(2) * t125 - t147) * t132 + (pkin(2) * t159 + t123) * t135) - g(3) * (pkin(2) * t160 + t121)) * MDP(7) + (-g(1) * (-t101 * t133 - t168) - g(2) * (t102 * t133 - t107) - t172) * MDP(20) + (-g(1) * (t101 * t130 - t167) - g(2) * (-t102 * t130 - t108) - t173) * MDP(21) + (-g(1) * (-t143 * t153 - t168) - g(2) * (t112 * t153 - t107) - t172) * MDP(22) + (-g(1) * (t143 * t155 - t167) - g(2) * (-t112 * t155 - t108) - t173) * MDP(23) + (-g(1) * (-pkin(5) * t168 - (t127 * t136 + t120) * t132 + (-t125 * t156 + t124) * t135 - t144 * t143) - g(2) * (-pkin(5) * t107 - (t125 * t136 - t147) * t132 + (t127 * t156 + t123) * t135 + t144 * t112) + (-t121 - ((t130 * pkin(5) + t136) * t135 + t144 * t132) * t126) * g(3)) * MDP(25) + (-MDP(3) + MDP(5)) * t90 + (-t131 * MDP(13) - t178 * t134 + MDP(4) - MDP(6)) * (-g(1) * t143 + g(2) * t112 + g(3) * t162); t148 * t90; t139 * MDP(13) + (-g(1) * (-t122 * t95 - t129 * t96) - g(2) * (t122 * t97 + t129 * t98) - g(3) * (-t115 * t122 - t116 * t129)) * MDP(25) + ((-t146 + t171) * MDP(20) + t139 * MDP(22)) * t133 + t178 * (g(1) * t96 - g(2) * t98 + g(3) * t116) + (-(MDP(21) + MDP(23)) * t171 + t146 * MDP(21) + t145 * MDP(23)) * t130; (-g(1) * (-t130 * t91 - t165) - g(2) * (t130 * t93 + t109) - t175) * MDP(20) + (-g(1) * (-t133 * t91 + t166) - g(2) * (t133 * t93 - t169) - t174) * MDP(21) + (-g(1) * (-t133 * t96 + t166) - g(2) * (t98 * t133 - t169) - t174) * MDP(23) + (pkin(5) * MDP(25) + MDP(22)) * (-g(1) * (-t130 * t96 - t165) - g(2) * (t130 * t98 + t109) - t175); -t139 * MDP(25);];
taug = t1;
