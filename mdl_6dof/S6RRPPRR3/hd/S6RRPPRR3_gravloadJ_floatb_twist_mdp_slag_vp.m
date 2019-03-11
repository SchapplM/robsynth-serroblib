% Calculate Gravitation load on the joints for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:39
% EndTime: 2019-03-09 08:58:41
% DurationCPUTime: 0.74s
% Computational Cost: add. (436->115), mult. (932->194), div. (0->0), fcn. (1148->14), ass. (0->57)
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t114 = pkin(12) + qJ(5);
t112 = sin(t114);
t113 = cos(t114);
t117 = sin(pkin(6));
t125 = cos(qJ(1));
t143 = t117 * t125;
t122 = sin(qJ(1));
t116 = sin(pkin(11));
t121 = sin(qJ(2));
t124 = cos(qJ(2));
t148 = cos(pkin(11));
t127 = -t121 * t116 + t124 * t148;
t105 = -t124 * t116 - t121 * t148;
t119 = cos(pkin(6));
t98 = t105 * t119;
t85 = -t122 * t127 + t125 * t98;
t81 = -t112 * t143 - t113 * t85;
t126 = t119 * t127;
t86 = t122 * t105 + t125 * t126;
t154 = t120 * t81 + t123 * t86;
t153 = -t120 * t86 + t123 * t81;
t139 = t122 * t124;
t141 = t121 * t125;
t102 = -t119 * t139 - t141;
t144 = t117 * t124;
t152 = -g(1) * t102 - g(3) * t144;
t147 = t113 * t120;
t146 = t113 * t123;
t145 = t117 * t122;
t140 = t122 * t121;
t138 = t124 * t125;
t136 = t119 * t138;
t111 = pkin(2) * t124 + pkin(1);
t99 = pkin(2) * t119 * t121 + (-pkin(8) - qJ(3)) * t117;
t134 = t125 * t111 - t122 * t99;
t132 = g(1) * t122 - g(2) * t125;
t88 = -t122 * t98 - t125 * t127;
t131 = -t122 * t111 - t125 * t99;
t130 = -t112 * t85 + t113 * t143;
t89 = t105 * t125 - t122 * t126;
t96 = t127 * t117;
t128 = g(1) * t89 + g(2) * t86 + g(3) * t96;
t118 = cos(pkin(12));
t115 = sin(pkin(12));
t107 = pkin(2) * t136;
t103 = -t119 * t140 + t138;
t101 = -t119 * t141 - t139;
t100 = -t136 + t140;
t97 = t105 * t117;
t92 = t112 * t119 - t113 * t97;
t84 = t112 * t145 - t113 * t88;
t83 = t112 * t88 + t113 * t145;
t79 = -t120 * t89 + t123 * t84;
t78 = -t120 * t84 - t123 * t89;
t1 = [t132 * MDP(2) + (-g(1) * t101 - g(2) * t103) * MDP(9) + (-g(1) * t100 - g(2) * t102) * MDP(10) + (-g(1) * t131 - g(2) * t134) * MDP(12) + (-g(1) * (t115 * t143 + t118 * t85) - g(2) * (t115 * t145 - t118 * t88)) * MDP(13) + (-g(1) * (-t115 * t85 + t118 * t143) - g(2) * (t115 * t88 + t118 * t145)) * MDP(14) + (-g(1) * t86 + g(2) * t89) * MDP(15) + (-g(1) * (pkin(3) * t85 + t86 * qJ(4) + t131) - g(2) * (-pkin(3) * t88 - qJ(4) * t89 + t134)) * MDP(16) + (g(1) * t81 - g(2) * t84) * MDP(22) + (-g(1) * t130 - g(2) * t83) * MDP(23) + (g(1) * t153 - g(2) * t79) * MDP(29) + (-g(1) * t154 - g(2) * t78) * MDP(30) + (-t117 * MDP(11) + MDP(3)) * (g(1) * t125 + g(2) * t122); (g(2) * t100 + t152) * MDP(9) + (g(3) * t117 * t121 + g(1) * t103 - g(2) * t101) * MDP(10) + (-g(2) * t107 + (g(2) * t140 + t152) * pkin(2)) * MDP(12) + (g(1) * t88 + g(2) * t85 + g(3) * t97) * MDP(15) + (-g(1) * (pkin(2) * t102 + t89 * pkin(3) - t88 * qJ(4)) - g(2) * (-pkin(2) * t140 + pkin(3) * t86 - qJ(4) * t85 + t107) - g(3) * (pkin(2) * t144 + pkin(3) * t96 - t97 * qJ(4))) * MDP(16) + (-g(1) * (-t120 * t88 + t146 * t89) - g(2) * (-t120 * t85 + t146 * t86) - g(3) * (-t120 * t97 + t146 * t96)) * MDP(29) + (-g(1) * (-t123 * t88 - t147 * t89) - g(2) * (-t123 * t85 - t147 * t86) - g(3) * (-t123 * t97 - t147 * t96)) * MDP(30) + (-MDP(13) * t118 + MDP(14) * t115 - MDP(22) * t113 + MDP(23) * t112) * t128; (MDP(12) + MDP(16)) * (-g(3) * t119 - t117 * t132); t128 * MDP(16); (g(1) * t84 + g(2) * t81 + g(3) * t92) * MDP(23) + (-MDP(29) * t123 + MDP(30) * t120 - MDP(22)) * (g(1) * t83 - g(2) * t130 + g(3) * (t112 * t97 + t113 * t119)); (-g(1) * t78 + g(2) * t154 - g(3) * (-t120 * t92 - t123 * t96)) * MDP(29) + (g(1) * t79 + g(2) * t153 - g(3) * (t120 * t96 - t123 * t92)) * MDP(30);];
taug  = t1;
