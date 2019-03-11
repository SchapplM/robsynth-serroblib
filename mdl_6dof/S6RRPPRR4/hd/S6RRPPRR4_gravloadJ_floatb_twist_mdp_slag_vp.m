% Calculate Gravitation load on the joints for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:04
% EndTime: 2019-03-09 09:05:06
% DurationCPUTime: 0.61s
% Computational Cost: add. (352->108), mult. (876->181), div. (0->0), fcn. (1076->12), ass. (0->54)
t120 = sin(pkin(11));
t125 = sin(qJ(2));
t129 = cos(qJ(2));
t158 = cos(pkin(11));
t112 = -t120 * t129 - t125 * t158;
t133 = -t120 * t125 + t129 * t158;
t122 = cos(pkin(6));
t130 = cos(qJ(1));
t146 = t130 * t125;
t126 = sin(qJ(1));
t149 = t126 * t129;
t108 = -t122 * t149 - t146;
t121 = sin(pkin(6));
t156 = t121 * t129;
t164 = -g(1) * t108 - g(3) * t156;
t123 = sin(qJ(6));
t127 = cos(qJ(6));
t124 = sin(qJ(5));
t128 = cos(qJ(5));
t155 = t121 * t130;
t131 = t133 * t122;
t93 = t126 * t112 + t130 * t131;
t135 = t124 * t93 + t128 * t155;
t144 = t112 * t122;
t92 = -t126 * t133 + t130 * t144;
t163 = t123 * t135 - t127 * t92;
t97 = t126 * t144 + t130 * t133;
t162 = t123 * t92 + t127 * t135;
t157 = t121 * t126;
t154 = t123 * t124;
t153 = t124 * t127;
t150 = t126 * t125;
t145 = t130 * t129;
t142 = t122 * t145;
t105 = t122 * t125 * pkin(2) + (-pkin(8) - qJ(3)) * t121;
t119 = pkin(2) * t129 + pkin(1);
t139 = -t105 * t126 + t119 * t130;
t137 = g(1) * t126 - g(2) * t130;
t136 = -t105 * t130 - t119 * t126;
t88 = t124 * t155 - t128 * t93;
t113 = pkin(2) * t142;
t109 = -t122 * t150 + t145;
t107 = -t122 * t146 - t149;
t106 = -t142 + t150;
t104 = t112 * t121;
t103 = t133 * t121;
t99 = -t103 * t124 + t122 * t128;
t96 = t112 * t130 - t126 * t131;
t87 = -t124 * t96 + t128 * t157;
t86 = -t124 * t157 - t128 * t96;
t85 = t123 * t97 + t127 * t87;
t84 = -t123 * t87 + t127 * t97;
t83 = g(1) * t96 + g(2) * t93 + g(3) * t103;
t1 = [t137 * MDP(2) + (-g(1) * t107 - g(2) * t109) * MDP(9) + (-g(1) * t106 - g(2) * t108) * MDP(10) + (-g(1) * t136 - g(2) * t139) * MDP(12) + (g(1) * t92 + g(2) * t97) * MDP(14) + (-g(1) * t93 + g(2) * t96) * MDP(15) + (-g(1) * (pkin(3) * t92 + qJ(4) * t93 + t136) - g(2) * (pkin(3) * t97 - qJ(4) * t96 + t139)) * MDP(16) + (-g(1) * t135 - g(2) * t87) * MDP(22) + (g(1) * t88 - g(2) * t86) * MDP(23) + (-g(1) * t162 - g(2) * t85) * MDP(29) + (g(1) * t163 - g(2) * t84) * MDP(30) + (MDP(3) - (MDP(11) + MDP(13)) * t121) * (g(1) * t130 + g(2) * t126); (g(2) * t106 + t164) * MDP(9) + (g(3) * t121 * t125 + g(1) * t109 - g(2) * t107) * MDP(10) + (-g(2) * t113 + (g(2) * t150 + t164) * pkin(2)) * MDP(12) + t83 * MDP(14) + (-g(1) * (pkin(2) * t108 + t96 * pkin(3) + qJ(4) * t97) - g(2) * (-pkin(2) * t150 + pkin(3) * t93 - qJ(4) * t92 + t113) - g(3) * (pkin(2) * t156 + pkin(3) * t103 - qJ(4) * t104)) * MDP(16) + (-g(1) * (t123 * t96 + t153 * t97) - g(2) * (t123 * t93 - t153 * t92) - g(3) * (t103 * t123 - t104 * t153)) * MDP(29) + (-g(1) * (t127 * t96 - t154 * t97) - g(2) * (t127 * t93 + t154 * t92) - g(3) * (t103 * t127 + t104 * t154)) * MDP(30) + (MDP(22) * t124 + MDP(23) * t128 + MDP(15)) * (-g(1) * t97 + g(2) * t92 + g(3) * t104); (MDP(12) + MDP(16)) * (-g(3) * t122 - t121 * t137); t83 * MDP(16); (g(1) * t87 - g(2) * t135 + g(3) * t99) * MDP(23) + (-MDP(29) * t127 + MDP(30) * t123 - MDP(22)) * (g(1) * t86 + g(2) * t88 + g(3) * (-t103 * t128 - t122 * t124)); (-g(1) * t84 - g(2) * t163 - g(3) * (-t104 * t127 - t123 * t99)) * MDP(29) + (g(1) * t85 - g(2) * t162 - g(3) * (t104 * t123 - t127 * t99)) * MDP(30);];
taug  = t1;
