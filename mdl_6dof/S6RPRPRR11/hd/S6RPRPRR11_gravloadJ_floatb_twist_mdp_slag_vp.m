% Calculate Gravitation load on the joints for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:15:14
% EndTime: 2019-03-09 04:15:18
% DurationCPUTime: 1.07s
% Computational Cost: add. (596->115), mult. (1460->195), div. (0->0), fcn. (1850->16), ass. (0->56)
t131 = sin(qJ(6));
t133 = cos(qJ(6));
t134 = cos(qJ(1));
t162 = sin(pkin(12));
t166 = cos(pkin(6));
t148 = t166 * t162;
t164 = cos(pkin(12));
t167 = sin(qJ(1));
t115 = t134 * t148 + t164 * t167;
t132 = sin(qJ(3));
t168 = cos(qJ(3));
t150 = t166 * t164;
t140 = -t134 * t150 + t167 * t162;
t129 = sin(pkin(6));
t163 = sin(pkin(7));
t156 = t129 * t163;
t165 = cos(pkin(7));
t175 = t134 * t156 + t140 * t165;
t102 = -t115 * t168 + t132 * t175;
t157 = t129 * t165;
t108 = -t134 * t157 + t140 * t163;
t127 = pkin(13) + qJ(5);
t124 = sin(t127);
t125 = cos(t127);
t94 = t102 * t125 - t108 * t124;
t99 = t115 * t132 + t168 * t175;
t179 = t131 * t94 + t133 * t99;
t178 = -t131 * t99 + t133 * t94;
t174 = t102 * t124 + t108 * t125;
t136 = t134 * t162 + t150 * t167;
t170 = t136 * t165 - t167 * t156;
t169 = MDP(14) - MDP(17);
t161 = qJ(2) * t129;
t160 = t125 * t131;
t159 = t125 * t133;
t158 = t134 * pkin(1) + t167 * t161;
t153 = -pkin(1) * t167 + t134 * t161;
t149 = t166 * t163;
t147 = t165 * t164;
t145 = g(1) * t167 - g(2) * t134;
t116 = t134 * t164 - t148 * t167;
t103 = t116 * t132 + t170 * t168;
t105 = -t149 * t168 + (t132 * t162 - t147 * t168) * t129;
t142 = g(1) * t103 + g(2) * t99 + g(3) * t105;
t109 = t136 * t163 + t157 * t167;
t130 = cos(pkin(13));
t128 = sin(pkin(13));
t114 = -t156 * t164 + t165 * t166;
t106 = t132 * t149 + (t132 * t147 + t168 * t162) * t129;
t104 = t116 * t168 - t170 * t132;
t98 = t106 * t125 + t114 * t124;
t96 = t104 * t125 + t109 * t124;
t95 = -t104 * t124 + t109 * t125;
t91 = t103 * t131 + t133 * t96;
t90 = t103 * t133 - t131 * t96;
t1 = [t145 * MDP(2) + (g(1) * t115 - g(2) * t116) * MDP(4) + (-g(1) * t140 + g(2) * t136) * MDP(5) + (-g(1) * t153 - g(2) * t158) * MDP(7) + (-g(1) * t102 - g(2) * t104) * MDP(13) + (-g(1) * (t102 * t130 - t108 * t128) - g(2) * (t104 * t130 + t109 * t128)) * MDP(15) + (-g(1) * (-t102 * t128 - t108 * t130) - g(2) * (-t104 * t128 + t109 * t130)) * MDP(16) + (-g(1) * (-t115 * pkin(2) + t102 * pkin(3) - qJ(4) * t99 + t153) - g(2) * (t116 * pkin(2) + t104 * pkin(3) + t103 * qJ(4) + t158) + (g(1) * t108 - g(2) * t109) * pkin(9)) * MDP(18) + (-g(1) * t94 - g(2) * t96) * MDP(24) + (g(1) * t174 - g(2) * t95) * MDP(25) + (-g(1) * t178 - g(2) * t91) * MDP(31) + (g(1) * t179 - g(2) * t90) * MDP(32) + t169 * (-g(1) * t99 + g(2) * t103) + (MDP(6) * t129 - MDP(3)) * (-g(1) * t134 - g(2) * t167); (MDP(18) + MDP(7)) * (-g(3) * t166 - t129 * t145); (-g(1) * (-pkin(3) * t103 + qJ(4) * t104) - g(2) * (-pkin(3) * t99 - qJ(4) * t102) - g(3) * (-pkin(3) * t105 + qJ(4) * t106)) * MDP(18) + (-g(1) * (-t103 * t159 + t104 * t131) - g(2) * (-t102 * t131 - t159 * t99) - g(3) * (-t105 * t159 + t106 * t131)) * MDP(31) + (-g(1) * (t103 * t160 + t104 * t133) - g(2) * (-t102 * t133 + t160 * t99) - g(3) * (t105 * t160 + t106 * t133)) * MDP(32) + t169 * (g(1) * t104 - g(2) * t102 + g(3) * t106) + (MDP(15) * t130 - MDP(16) * t128 + MDP(24) * t125 - MDP(25) * t124 + MDP(13)) * t142; -t142 * MDP(18); (g(1) * t96 - g(2) * t94 + g(3) * t98) * MDP(25) + (-MDP(31) * t133 + MDP(32) * t131 - MDP(24)) * (g(1) * t95 + g(2) * t174 + g(3) * (-t106 * t124 + t114 * t125)); (-g(1) * t90 - g(2) * t179 - g(3) * (t105 * t133 - t131 * t98)) * MDP(31) + (g(1) * t91 - g(2) * t178 - g(3) * (-t105 * t131 - t133 * t98)) * MDP(32);];
taug  = t1;
