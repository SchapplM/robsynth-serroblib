% Calculate Gravitation load on the joints for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:22
% EndTime: 2019-03-09 04:25:25
% DurationCPUTime: 0.77s
% Computational Cost: add. (503->108), mult. (1375->182), div. (0->0), fcn. (1740->14), ass. (0->59)
t125 = sin(qJ(6));
t129 = cos(qJ(6));
t124 = cos(pkin(6));
t119 = sin(pkin(12));
t128 = sin(qJ(1));
t149 = t128 * t119;
t122 = cos(pkin(12));
t131 = cos(qJ(1));
t153 = t122 * t131;
t109 = -t124 * t153 + t149;
t120 = sin(pkin(7));
t123 = cos(pkin(7));
t121 = sin(pkin(6));
t154 = t121 * t131;
t101 = -t109 * t120 + t123 * t154;
t126 = sin(qJ(5));
t130 = cos(qJ(5));
t147 = t131 * t119;
t148 = t128 * t122;
t110 = t124 * t147 + t148;
t127 = sin(qJ(3));
t158 = cos(qJ(3));
t145 = t120 * t158;
t142 = t121 * t145;
t144 = t123 * t158;
t93 = t109 * t144 + t110 * t127 + t131 * t142;
t86 = t101 * t130 - t126 * t93;
t152 = t123 * t127;
t157 = t120 * t127;
t94 = -t109 * t152 + t110 * t158 - t154 * t157;
t165 = t125 * t86 + t129 * t94;
t164 = -t125 * t94 + t129 * t86;
t163 = MDP(13) - MDP(16);
t160 = t101 * t126 + t130 * t93;
t159 = MDP(14) - MDP(17);
t156 = t121 * t122;
t155 = t121 * t128;
t151 = t125 * t126;
t150 = t126 * t129;
t146 = t131 * pkin(1) + qJ(2) * t155;
t143 = -t128 * pkin(1) + qJ(2) * t154;
t136 = t124 * t148 + t147;
t103 = t120 * t136 + t123 * t155;
t139 = -g(1) * t101 - g(2) * t103;
t137 = g(1) * t128 - g(2) * t131;
t111 = -t124 * t149 + t153;
t132 = t136 * t123;
t97 = t111 * t127 - t128 * t142 + t132 * t158;
t99 = t119 * t121 * t127 - t124 * t145 - t144 * t156;
t134 = g(1) * t97 + g(2) * t93 + g(3) * t99;
t108 = -t120 * t156 + t123 * t124;
t100 = t124 * t157 + (t119 * t158 + t122 * t152) * t121;
t98 = t111 * t158 + (t120 * t155 - t132) * t127;
t92 = t108 * t130 + t126 * t99;
t90 = t103 * t130 + t126 * t97;
t89 = -t103 * t126 + t130 * t97;
t84 = t125 * t98 + t129 * t90;
t83 = -t125 * t90 + t129 * t98;
t1 = [t137 * MDP(2) + (g(1) * t110 - g(2) * t111) * MDP(4) + (-g(1) * t109 + g(2) * t136) * MDP(5) + (-g(1) * t143 - g(2) * t146) * MDP(7) + t139 * MDP(15) + (-g(1) * (-t110 * pkin(2) - pkin(3) * t94 - qJ(4) * t93 + t143) - g(2) * (t111 * pkin(2) + t98 * pkin(3) + t97 * qJ(4) + t146) + t139 * pkin(9)) * MDP(18) + (-g(1) * t86 - g(2) * t90) * MDP(24) + (g(1) * t160 - g(2) * t89) * MDP(25) + (-g(1) * t164 - g(2) * t84) * MDP(31) + (g(1) * t165 - g(2) * t83) * MDP(32) + t159 * (-g(1) * t93 + g(2) * t97) - t163 * (-g(1) * t94 + g(2) * t98) + (-MDP(6) * t121 + MDP(3)) * (g(1) * t131 + g(2) * t128); (MDP(18) + MDP(7)) * (-g(3) * t124 - t121 * t137); (-g(1) * (-pkin(3) * t97 + qJ(4) * t98) - g(2) * (-pkin(3) * t93 + qJ(4) * t94) - g(3) * (-pkin(3) * t99 + qJ(4) * t100)) * MDP(18) + (-g(1) * (-t125 * t97 + t150 * t98) - g(2) * (-t125 * t93 + t150 * t94) - g(3) * (t100 * t150 - t125 * t99)) * MDP(31) + (-g(1) * (-t129 * t97 - t151 * t98) - g(2) * (-t93 * t129 - t151 * t94) - g(3) * (-t100 * t151 - t129 * t99)) * MDP(32) + (-MDP(24) * t126 - MDP(25) * t130 + t159) * (g(1) * t98 + g(2) * t94 + g(3) * t100) + t163 * t134; -t134 * MDP(18); (g(1) * t90 - g(2) * t86 + g(3) * t92) * MDP(25) + (-MDP(31) * t129 + MDP(32) * t125 - MDP(24)) * (g(1) * t89 + g(2) * t160 + g(3) * (-t108 * t126 + t130 * t99)); (-g(1) * t83 - g(2) * t165 - g(3) * (t100 * t129 - t125 * t92)) * MDP(31) + (g(1) * t84 - g(2) * t164 - g(3) * (-t100 * t125 - t129 * t92)) * MDP(32);];
taug  = t1;
