% Calculate Gravitation load on the joints for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:34
% EndTime: 2019-03-09 13:01:36
% DurationCPUTime: 0.70s
% Computational Cost: add. (324->114), mult. (771->176), div. (0->0), fcn. (889->10), ass. (0->53)
t165 = MDP(9) - MDP(12);
t164 = -MDP(10) + MDP(13);
t161 = MDP(21) - MDP(29);
t118 = sin(qJ(2));
t119 = sin(qJ(1));
t122 = cos(qJ(2));
t123 = cos(qJ(1));
t154 = cos(pkin(6));
t138 = t123 * t154;
t102 = t118 * t138 + t119 * t122;
t139 = t119 * t154;
t104 = -t118 * t139 + t123 * t122;
t162 = -g(1) * t104 - g(2) * t102;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t101 = t119 * t118 - t122 * t138;
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t114 = sin(pkin(6));
t146 = t114 * t123;
t131 = -t101 * t117 + t121 * t146;
t160 = t102 * t120 + t116 * t131;
t159 = -t102 * t116 + t120 * t131;
t155 = g(3) * t114;
t153 = t101 * t116;
t103 = t123 * t118 + t122 * t139;
t150 = t103 * t116;
t149 = t114 * t118;
t148 = t114 * t119;
t147 = t114 * t122;
t145 = t116 * t117;
t144 = t116 * t118;
t143 = t117 * t120;
t142 = t118 * t120;
t141 = pkin(5) * t116 + pkin(9);
t140 = g(3) * (pkin(2) * t147 + qJ(3) * t149);
t87 = t103 * t117 + t121 * t148;
t82 = t104 * t120 - t87 * t116;
t112 = t120 * pkin(5) + pkin(4);
t115 = -qJ(6) - pkin(10);
t133 = t112 * t117 + t115 * t121;
t132 = t123 * pkin(1) + t104 * pkin(2) + pkin(8) * t148 + t103 * qJ(3);
t130 = t101 * t121 + t117 * t146;
t86 = -t103 * t121 + t117 * t148;
t99 = t117 * t154 + t121 * t147;
t128 = g(1) * t86 - g(2) * t130 + g(3) * t99;
t126 = -t119 * pkin(1) - t102 * pkin(2) + pkin(8) * t146 - t101 * qJ(3);
t85 = -g(1) * t103 - g(2) * t101 + g(3) * t147;
t100 = -t117 * t147 + t121 * t154;
t97 = t103 * pkin(2);
t95 = t101 * pkin(2);
t83 = t104 * t116 + t87 * t120;
t1 = [(g(1) * t119 - g(2) * t123) * MDP(2) + (-g(1) * t126 - g(2) * t132) * MDP(14) + (-g(1) * t131 - g(2) * t87) * MDP(20) + (-g(1) * t159 - g(2) * t83) * MDP(27) + (g(1) * t160 - g(2) * t82) * MDP(28) + (-g(1) * (pkin(3) * t146 - t102 * t141 + t112 * t131 - t115 * t130 + t126) - g(2) * (pkin(3) * t148 + t104 * t141 + t87 * t112 - t86 * t115 + t132)) * MDP(30) + t161 * (g(1) * t130 + g(2) * t86) + t164 * (g(1) * t101 - g(2) * t103) + t165 * (g(1) * t102 - g(2) * t104) + (-t114 * MDP(11) + MDP(3)) * (g(1) * t123 + g(2) * t119); (-g(1) * (t104 * qJ(3) - t97) - g(2) * (t102 * qJ(3) - t95) - t140) * MDP(14) + (-g(1) * (t104 * t143 - t150) - g(2) * (t102 * t143 - t153) - (t116 * t122 + t117 * t142) * t155) * MDP(27) + (-g(1) * (-t103 * t120 - t104 * t145) - g(2) * (-t101 * t120 - t102 * t145) - (-t117 * t144 + t120 * t122) * t155) * MDP(28) + (-g(1) * (-pkin(5) * t150 - t103 * pkin(9) - t97) - g(2) * (-pkin(5) * t153 - t101 * pkin(9) - t95) - t140 - (t118 * t133 + t122 * t141) * t155 + t162 * (qJ(3) + t133)) * MDP(30) - t165 * t85 + (-t117 * MDP(20) - t161 * t121 - t164) * (g(3) * t149 - t162); (MDP(14) + MDP(30)) * t85; (-g(1) * (-t86 * t112 - t87 * t115) - g(2) * (t112 * t130 + t115 * t131) - g(3) * (-t100 * t115 - t99 * t112)) * MDP(30) + t161 * (g(1) * t87 - g(2) * t131 + g(3) * t100) + (MDP(27) * t120 - MDP(28) * t116 + MDP(20)) * t128; (g(1) * t83 - g(2) * t159 - g(3) * (-t100 * t120 - t114 * t144)) * MDP(28) + (pkin(5) * MDP(30) + MDP(27)) * (-g(2) * t160 - g(3) * (-t100 * t116 + t114 * t142) - g(1) * t82); -t128 * MDP(30);];
taug  = t1;
