% Calculate Gravitation load on the joints for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:46:07
% EndTime: 2019-03-09 07:46:10
% DurationCPUTime: 1.10s
% Computational Cost: add. (707->116), mult. (1820->206), div. (0->0), fcn. (2353->16), ass. (0->59)
t134 = cos(qJ(1));
t160 = sin(pkin(13));
t164 = cos(pkin(6));
t148 = t164 * t160;
t162 = cos(pkin(13));
t166 = sin(qJ(1));
t118 = t134 * t148 + t166 * t162;
t131 = sin(qJ(3));
t167 = cos(qJ(3));
t149 = t164 * t162;
t140 = -t134 * t149 + t166 * t160;
t128 = sin(pkin(6));
t161 = sin(pkin(7));
t153 = t128 * t161;
t163 = cos(pkin(7));
t173 = t134 * t153 + t140 * t163;
t104 = t118 * t131 + t173 * t167;
t127 = qJ(5) + qJ(6);
t125 = sin(t127);
t126 = cos(t127);
t107 = -t118 * t167 + t173 * t131;
t154 = t128 * t163;
t113 = -t134 * t154 + t140 * t161;
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t99 = t107 * t133 - t113 * t130;
t181 = t104 * t126 + t99 * t125;
t180 = -t104 * t125 + t99 * t126;
t129 = sin(qJ(5));
t132 = cos(qJ(5));
t179 = t104 * t132 + t99 * t129;
t178 = -t104 * t129 + t99 * t132;
t172 = t107 * t130 + t113 * t133;
t136 = t134 * t160 + t166 * t149;
t168 = t136 * t163 - t166 * t153;
t146 = t161 * t164;
t147 = t163 * t162;
t111 = t131 * t146 + (t131 * t147 + t167 * t160) * t128;
t117 = -t162 * t153 + t164 * t163;
t103 = t111 * t133 + t117 * t130;
t110 = -t167 * t146 + (t131 * t160 - t147 * t167) * t128;
t119 = t134 * t162 - t166 * t148;
t109 = t119 * t167 - t168 * t131;
t114 = t136 * t161 + t166 * t154;
t101 = t109 * t133 + t114 * t130;
t108 = t119 * t131 + t168 * t167;
t93 = -t101 * t125 + t108 * t126;
t94 = t101 * t126 + t108 * t125;
t165 = (-g(1) * t93 - g(2) * t181 - g(3) * (-t103 * t125 + t110 * t126)) * MDP(34) + (g(1) * t94 - g(2) * t180 - g(3) * (-t103 * t126 - t110 * t125)) * MDP(35);
t159 = qJ(2) * t128;
t158 = t125 * t133;
t157 = t126 * t133;
t156 = t129 * t133;
t155 = t132 * t133;
t144 = g(1) * t166 - g(2) * t134;
t100 = -t109 * t130 + t114 * t133;
t96 = t101 * t132 + t108 * t129;
t95 = -t101 * t129 + t108 * t132;
t1 = [t144 * MDP(2) + (g(1) * t118 - g(2) * t119) * MDP(4) + (-g(1) * t140 + g(2) * t136) * MDP(5) + (-g(1) * (-t166 * pkin(1) + t134 * t159) - g(2) * (t134 * pkin(1) + t166 * t159)) * MDP(7) + (-g(1) * t107 - g(2) * t109) * MDP(13) + (-g(1) * t104 + g(2) * t108) * MDP(14) + (-g(1) * t99 - g(2) * t101) * MDP(20) + (g(1) * t172 - g(2) * t100) * MDP(21) + (-g(1) * t178 - g(2) * t96) * MDP(27) + (g(1) * t179 - g(2) * t95) * MDP(28) + (-g(1) * t180 - g(2) * t94) * MDP(34) + (g(1) * t181 - g(2) * t93) * MDP(35) + (t128 * MDP(6) - MDP(3)) * (-g(1) * t134 - g(2) * t166); (-g(3) * t164 - t144 * t128) * MDP(7); (g(1) * t109 - g(2) * t107 + g(3) * t111) * MDP(14) + (-g(1) * (-t108 * t155 + t109 * t129) - g(2) * (-t104 * t155 - t107 * t129) - g(3) * (-t110 * t155 + t111 * t129)) * MDP(27) + (-g(1) * (t108 * t156 + t109 * t132) - g(2) * (t104 * t156 - t107 * t132) - g(3) * (t110 * t156 + t111 * t132)) * MDP(28) + (-g(1) * (-t108 * t157 + t109 * t125) - g(2) * (-t104 * t157 - t107 * t125) - g(3) * (-t110 * t157 + t111 * t125)) * MDP(34) + (-g(1) * (t108 * t158 + t109 * t126) - g(2) * (t104 * t158 - t107 * t126) - g(3) * (t110 * t158 + t111 * t126)) * MDP(35) + (t133 * MDP(20) - MDP(21) * t130 + MDP(13)) * (g(1) * t108 + g(2) * t104 + g(3) * t110); (g(1) * t101 - g(2) * t99 + g(3) * t103) * MDP(21) + (-MDP(27) * t132 + MDP(28) * t129 - MDP(34) * t126 + MDP(35) * t125 - MDP(20)) * (g(1) * t100 + g(2) * t172 + g(3) * (-t111 * t130 + t117 * t133)); (-g(1) * t95 - g(2) * t179 - g(3) * (-t103 * t129 + t110 * t132)) * MDP(27) + (g(1) * t96 - g(2) * t178 - g(3) * (-t103 * t132 - t110 * t129)) * MDP(28) + t165; t165;];
taug  = t1;
