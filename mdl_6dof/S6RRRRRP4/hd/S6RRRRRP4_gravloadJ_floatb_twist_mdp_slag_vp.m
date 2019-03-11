% Calculate Gravitation load on the joints for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:47
% EndTime: 2019-03-10 01:15:48
% DurationCPUTime: 0.41s
% Computational Cost: add. (514->85), mult. (532->118), div. (0->0), fcn. (515->10), ass. (0->51)
t122 = qJ(4) + qJ(5);
t117 = sin(t122);
t119 = cos(t122);
t171 = pkin(5) * t119 + qJ(6) * t117;
t170 = MDP(17) - MDP(33);
t169 = MDP(30) + MDP(32);
t168 = MDP(31) - MDP(34);
t123 = qJ(2) + qJ(3);
t118 = sin(t123);
t128 = cos(qJ(2));
t130 = -pkin(10) - pkin(9);
t167 = t128 * pkin(2) - t118 * t130;
t126 = sin(qJ(1));
t129 = cos(qJ(1));
t141 = g(1) * t129 + g(2) * t126;
t120 = cos(t123);
t127 = cos(qJ(4));
t148 = t127 * t129;
t124 = sin(qJ(4));
t152 = t124 * t126;
t105 = t120 * t152 + t148;
t149 = t126 * t127;
t151 = t124 * t129;
t107 = -t120 * t151 + t149;
t160 = g(3) * t118;
t166 = -g(1) * t107 + g(2) * t105 + t124 * t160;
t157 = t117 * t118;
t156 = t117 * t126;
t155 = t118 * t119;
t153 = t119 * t129;
t115 = pkin(4) * t127 + pkin(3);
t111 = t120 * t115;
t150 = t126 * t119;
t147 = t129 * t117;
t145 = t171 * t120 + t111;
t144 = pkin(4) * t124 + pkin(7) + pkin(8);
t101 = t120 * t150 - t147;
t103 = t120 * t153 + t156;
t100 = t120 * t156 + t153;
t102 = t120 * t147 - t150;
t84 = g(1) * t102 + g(2) * t100 + g(3) * t157;
t143 = t169 * t84 + t168 * (g(1) * t103 + g(2) * t101 + g(3) * t155);
t139 = t115 + t171;
t138 = t141 * t120;
t137 = t111 + pkin(1) + t167;
t136 = t170 * (t138 + t160) + (MDP(23) * t127 - MDP(24) * t124 - t168 * t117 + t169 * t119 + MDP(16)) * (-g(3) * t120 + t141 * t118);
t132 = -g(1) * (-t102 * pkin(5) + qJ(6) * t103) - g(2) * (-t100 * pkin(5) + qJ(6) * t101) - g(3) * (-pkin(5) * t157 + qJ(6) * t155);
t125 = sin(qJ(2));
t108 = t120 * t148 + t152;
t106 = -t120 * t149 + t151;
t1 = [t141 * MDP(3) + (-g(1) * t106 - g(2) * t108) * MDP(23) + (-g(1) * t105 - g(2) * t107) * MDP(24) + (-g(1) * (-t101 * pkin(5) - t100 * qJ(6)) - g(2) * (t103 * pkin(5) + t102 * qJ(6)) + (-g(1) * t144 - g(2) * t137) * t129 + (g(1) * t137 - g(2) * t144) * t126) * MDP(35) + t169 * (g(1) * t101 - g(2) * t103) - t168 * (g(1) * t100 - g(2) * t102) + (-t125 * MDP(10) + t120 * MDP(16) + t128 * MDP(9) - t170 * t118 + MDP(2)) * (g(1) * t126 - g(2) * t129); (-g(3) * t128 + t141 * t125) * MDP(9) + (g(3) * t125 + t141 * t128) * MDP(10) + (-g(3) * (t145 + t167) + t141 * (pkin(2) * t125 + t139 * t118 + t120 * t130)) * MDP(35) + t136; (-g(3) * t145 + t130 * t138 + (g(3) * t130 + t139 * t141) * t118) * MDP(35) + t136; t166 * MDP(23) + (g(1) * t108 - g(2) * t106 + t127 * t160) * MDP(24) + (pkin(4) * t166 + t132) * MDP(35) + t143; t132 * MDP(35) + t143; -t84 * MDP(35);];
taug  = t1;
