% Calculate Gravitation load on the joints for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:14
% EndTime: 2019-12-31 22:55:17
% DurationCPUTime: 0.79s
% Computational Cost: add. (477->122), mult. (1324->222), div. (0->0), fcn. (1696->14), ass. (0->66)
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t132 = sin(qJ(2));
t135 = cos(qJ(2));
t136 = cos(qJ(1));
t159 = cos(pkin(5));
t142 = t136 * t159;
t160 = sin(qJ(1));
t116 = t160 * t132 - t135 * t142;
t126 = sin(pkin(6));
t128 = cos(pkin(6));
t127 = sin(pkin(5));
t154 = t127 * t136;
t108 = -t116 * t126 + t128 * t154;
t130 = sin(qJ(4));
t134 = cos(qJ(4));
t117 = t132 * t142 + t160 * t135;
t131 = sin(qJ(3));
t147 = t126 * t154;
t153 = t128 * t131;
t161 = cos(qJ(3));
t96 = -t116 * t153 + t117 * t161 - t131 * t147;
t88 = t108 * t130 - t96 * t134;
t145 = t128 * t161;
t95 = t116 * t145 + t117 * t131 + t161 * t147;
t166 = t88 * t129 + t95 * t133;
t165 = -t95 * t129 + t88 * t133;
t162 = t108 * t134 + t96 * t130;
t158 = t126 * t130;
t157 = t126 * t134;
t156 = t127 * t132;
t155 = t127 * t135;
t152 = t129 * t134;
t151 = t131 * t132;
t150 = t131 * t135;
t149 = t133 * t134;
t148 = t126 * t156;
t146 = t127 * t160;
t144 = t161 * t132;
t143 = t161 * t135;
t141 = t159 * t126;
t140 = t126 * t146;
t139 = t159 * t160;
t119 = -t132 * t139 + t136 * t135;
t118 = -t136 * t132 - t135 * t139;
t115 = -t126 * t155 + t159 * t128;
t114 = (-t128 * t151 + t143) * t127;
t113 = (t128 * t144 + t150) * t127;
t110 = -t118 * t126 + t128 * t146;
t107 = t131 * t141 + (t128 * t150 + t144) * t127;
t106 = -t161 * t141 + (-t128 * t143 + t151) * t127;
t105 = t114 * t134 + t130 * t148;
t104 = t118 * t161 - t119 * t153;
t103 = t118 * t131 + t119 * t145;
t102 = -t116 * t161 - t117 * t153;
t101 = -t116 * t131 + t117 * t145;
t100 = t119 * t161 + (t118 * t128 + t140) * t131;
t99 = -t118 * t145 + t119 * t131 - t161 * t140;
t94 = t107 * t134 + t115 * t130;
t92 = t104 * t134 + t119 * t158;
t91 = t102 * t134 + t117 * t158;
t90 = t100 * t134 + t110 * t130;
t89 = -t100 * t130 + t110 * t134;
t85 = t99 * t129 + t90 * t133;
t84 = -t90 * t129 + t99 * t133;
t1 = [(g(1) * t160 - g(2) * t136) * MDP(2) + (g(1) * t136 + g(2) * t160) * MDP(3) + (g(1) * t117 - g(2) * t119) * MDP(9) + (-g(1) * t116 - g(2) * t118) * MDP(10) + (g(1) * t96 - g(2) * t100) * MDP(16) + (-g(1) * t95 + g(2) * t99) * MDP(17) + (-g(1) * t88 - g(2) * t90) * MDP(23) + (-g(1) * t162 - g(2) * t89) * MDP(24) + (-g(1) * t165 - g(2) * t85) * MDP(30) + (g(1) * t166 - g(2) * t84) * MDP(31); (-g(1) * t118 + g(2) * t116 - g(3) * t155) * MDP(9) + (g(1) * t119 + g(2) * t117 + g(3) * t156) * MDP(10) + (-g(1) * t104 - g(2) * t102 - g(3) * t114) * MDP(16) + (g(1) * t103 + g(2) * t101 + g(3) * t113) * MDP(17) + (-g(1) * t92 - g(2) * t91 - g(3) * t105) * MDP(23) + (-g(1) * (-t104 * t130 + t119 * t157) - g(2) * (-t102 * t130 + t117 * t157) - g(3) * (-t114 * t130 + t134 * t148)) * MDP(24) + (-g(1) * (t103 * t129 + t92 * t133) - g(2) * (t101 * t129 + t91 * t133) - g(3) * (t105 * t133 + t113 * t129)) * MDP(30) + (-g(1) * (t103 * t133 - t92 * t129) - g(2) * (t101 * t133 - t91 * t129) - g(3) * (-t105 * t129 + t113 * t133)) * MDP(31); (g(1) * t100 + g(2) * t96 + g(3) * t107) * MDP(17) + (-g(1) * (t100 * t129 - t99 * t149) - g(2) * (t96 * t129 - t95 * t149) - g(3) * (-t106 * t149 + t107 * t129)) * MDP(30) + (-g(1) * (t100 * t133 + t99 * t152) - g(2) * (t96 * t133 + t95 * t152) - g(3) * (t106 * t152 + t107 * t133)) * MDP(31) + (t134 * MDP(23) - MDP(24) * t130 + MDP(16)) * (g(1) * t99 + g(2) * t95 + g(3) * t106); (g(1) * t90 - g(2) * t88 + g(3) * t94) * MDP(24) + (-MDP(30) * t133 + MDP(31) * t129 - MDP(23)) * (g(1) * t89 - g(2) * t162 + g(3) * (-t107 * t130 + t115 * t134)); (-g(1) * t84 - g(2) * t166 - g(3) * (t106 * t133 - t94 * t129)) * MDP(30) + (g(1) * t85 - g(2) * t165 - g(3) * (-t106 * t129 - t94 * t133)) * MDP(31);];
taug = t1;
