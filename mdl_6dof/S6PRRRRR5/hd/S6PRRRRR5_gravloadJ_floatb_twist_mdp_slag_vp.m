% Calculate Gravitation load on the joints for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:03
% EndTime: 2019-03-09 01:10:05
% DurationCPUTime: 0.66s
% Computational Cost: add. (696->141), mult. (1803->254), div. (0->0), fcn. (2328->16), ass. (0->62)
t115 = sin(pkin(7));
t121 = sin(qJ(2));
t124 = cos(qJ(2));
t117 = cos(pkin(13));
t151 = cos(pkin(6));
t135 = t117 * t151;
t149 = sin(pkin(13));
t127 = t121 * t149 - t124 * t135;
t116 = sin(pkin(6));
t144 = t116 * t117;
t150 = cos(pkin(7));
t155 = t115 * t144 + t127 * t150;
t131 = t151 * t149;
t128 = t117 * t121 + t124 * t131;
t136 = t116 * t149;
t154 = -t115 * t136 + t128 * t150;
t153 = cos(qJ(3));
t114 = qJ(5) + qJ(6);
t112 = sin(t114);
t113 = cos(t114);
t119 = sin(qJ(4));
t123 = cos(qJ(4));
t106 = t121 * t135 + t124 * t149;
t120 = sin(qJ(3));
t87 = t106 * t153 - t155 * t120;
t99 = t115 * t127 - t144 * t150;
t81 = t119 * t99 + t123 * t87;
t100 = t115 * t128 + t136 * t150;
t107 = t117 * t124 - t121 * t131;
t89 = t107 * t153 - t154 * t120;
t83 = t100 * t119 + t123 * t89;
t86 = t106 * t120 + t155 * t153;
t88 = t107 * t120 + t154 * t153;
t142 = t116 * t124;
t105 = -t115 * t142 + t150 * t151;
t134 = t120 * t150;
t137 = t115 * t151;
t98 = t120 * t137 + (t121 * t153 + t124 * t134) * t116;
t91 = t105 * t119 + t123 * t98;
t132 = t150 * t153;
t143 = t116 * t121;
t97 = t120 * t143 - t132 * t142 - t137 * t153;
t152 = (-g(1) * (-t112 * t83 + t113 * t88) - g(2) * (-t112 * t81 + t113 * t86) - g(3) * (-t112 * t91 + t113 * t97)) * MDP(31) + (-g(1) * (-t112 * t88 - t113 * t83) - g(2) * (-t112 * t86 - t113 * t81) - g(3) * (-t112 * t97 - t113 * t91)) * MDP(32);
t148 = t112 * t123;
t147 = t113 * t123;
t146 = t115 * t119;
t145 = t115 * t123;
t118 = sin(qJ(5));
t141 = t118 * t123;
t122 = cos(qJ(5));
t140 = t122 * t123;
t139 = t115 * t143;
t104 = (-t121 * t134 + t124 * t153) * t116;
t103 = (t120 * t124 + t121 * t132) * t116;
t96 = t104 * t123 + t119 * t139;
t95 = -t107 * t134 - t128 * t153;
t94 = t107 * t132 - t120 * t128;
t93 = -t106 * t134 - t127 * t153;
t92 = t106 * t132 - t120 * t127;
t85 = t107 * t146 + t123 * t95;
t84 = t106 * t146 + t123 * t93;
t1 = [-g(3) * MDP(1); (g(1) * t128 + g(2) * t127 - g(3) * t142) * MDP(3) + (g(1) * t107 + g(2) * t106 + g(3) * t143) * MDP(4) + (-g(1) * t95 - g(2) * t93 - g(3) * t104) * MDP(10) + (g(1) * t94 + g(2) * t92 + g(3) * t103) * MDP(11) + (-g(1) * t85 - g(2) * t84 - g(3) * t96) * MDP(17) + (-g(1) * (t107 * t145 - t119 * t95) - g(2) * (t106 * t145 - t119 * t93) - g(3) * (-t104 * t119 + t123 * t139)) * MDP(18) + (-g(1) * (t118 * t94 + t122 * t85) - g(2) * (t118 * t92 + t122 * t84) - g(3) * (t103 * t118 + t122 * t96)) * MDP(24) + (-g(1) * (-t118 * t85 + t122 * t94) - g(2) * (-t118 * t84 + t122 * t92) - g(3) * (t103 * t122 - t118 * t96)) * MDP(25) + (-g(1) * (t112 * t94 + t113 * t85) - g(2) * (t112 * t92 + t113 * t84) - g(3) * (t103 * t112 + t113 * t96)) * MDP(31) + (-g(1) * (-t112 * t85 + t113 * t94) - g(2) * (-t112 * t84 + t113 * t92) - g(3) * (t103 * t113 - t112 * t96)) * MDP(32); (g(1) * t89 + g(2) * t87 + g(3) * t98) * MDP(11) + (-g(1) * (t118 * t89 - t140 * t88) - g(2) * (t118 * t87 - t140 * t86) - g(3) * (t118 * t98 - t140 * t97)) * MDP(24) + (-g(1) * (t122 * t89 + t141 * t88) - g(2) * (t122 * t87 + t141 * t86) - g(3) * (t122 * t98 + t141 * t97)) * MDP(25) + (-g(1) * (t112 * t89 - t147 * t88) - g(2) * (t112 * t87 - t147 * t86) - g(3) * (t112 * t98 - t147 * t97)) * MDP(31) + (-g(1) * (t113 * t89 + t148 * t88) - g(2) * (t113 * t87 + t148 * t86) - g(3) * (t113 * t98 + t148 * t97)) * MDP(32) + (MDP(17) * t123 - MDP(18) * t119 + MDP(10)) * (g(1) * t88 + g(2) * t86 + g(3) * t97); (g(1) * t83 + g(2) * t81 + g(3) * t91) * MDP(18) + (-MDP(24) * t122 + MDP(25) * t118 - MDP(31) * t113 + MDP(32) * t112 - MDP(17)) * (g(1) * (t100 * t123 - t119 * t89) + g(2) * (-t119 * t87 + t123 * t99) + g(3) * (t105 * t123 - t119 * t98)); (-g(1) * (-t118 * t83 + t122 * t88) - g(2) * (-t118 * t81 + t122 * t86) - g(3) * (-t118 * t91 + t122 * t97)) * MDP(24) + (-g(1) * (-t118 * t88 - t122 * t83) - g(2) * (-t118 * t86 - t122 * t81) - g(3) * (-t118 * t97 - t122 * t91)) * MDP(25) + t152; t152;];
taug  = t1;
