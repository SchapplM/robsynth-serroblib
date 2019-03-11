% Calculate Gravitation load on the joints for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:32
% EndTime: 2019-03-09 16:41:34
% DurationCPUTime: 0.59s
% Computational Cost: add. (483->100), mult. (509->135), div. (0->0), fcn. (470->10), ass. (0->53)
t116 = pkin(10) + qJ(5);
t111 = sin(t116);
t112 = cos(t116);
t161 = pkin(5) * t112 + qJ(6) * t111;
t158 = MDP(27) + MDP(29);
t157 = MDP(28) - MDP(31);
t160 = MDP(17) - MDP(20) - MDP(30);
t117 = qJ(2) + qJ(3);
t113 = sin(t117);
t114 = cos(t117);
t159 = t114 * pkin(3) + t113 * qJ(4);
t122 = sin(qJ(1));
t124 = cos(qJ(1));
t134 = g(1) * t124 + g(2) * t122;
t91 = -g(3) * t114 + t134 * t113;
t121 = sin(qJ(2));
t156 = pkin(2) * t121;
t155 = pkin(3) * t113;
t151 = g(3) * t113;
t120 = -pkin(9) - qJ(4);
t148 = t113 * t120;
t119 = cos(pkin(10));
t107 = pkin(4) * t119 + pkin(3);
t99 = t114 * t107;
t147 = t114 * t122;
t146 = t114 * t124;
t118 = sin(pkin(10));
t145 = t118 * t122;
t144 = t118 * t124;
t143 = t119 * t122;
t142 = t119 * t124;
t141 = t122 * t112;
t140 = t124 * t111;
t138 = t114 * t161 + t99;
t125 = -pkin(8) - pkin(7);
t137 = pkin(4) * t118 - t125;
t135 = -t155 - t156;
t131 = t99 - t148;
t130 = t107 + t161;
t129 = t134 * t114;
t93 = t111 * t147 + t112 * t124;
t95 = t114 * t140 - t141;
t77 = g(1) * t95 + g(2) * t93 + t111 * t151;
t127 = t160 * (t129 + t151) + (MDP(18) * t119 - MDP(19) * t118 - t157 * t111 + t158 * t112 + MDP(16)) * t91;
t123 = cos(qJ(2));
t115 = t123 * pkin(2);
t110 = t115 + pkin(1);
t103 = t124 * t110;
t102 = qJ(4) * t146;
t101 = qJ(4) * t147;
t96 = t111 * t122 + t112 * t146;
t94 = t114 * t141 - t140;
t1 = [t134 * MDP(3) + (-g(1) * (-t114 * t143 + t144) - g(2) * (t114 * t142 + t145)) * MDP(18) + (-g(1) * (t114 * t145 + t142) - g(2) * (-t114 * t144 + t143)) * MDP(19) + (-g(2) * t103 + (g(1) * t125 - g(2) * t159) * t124 + (-g(1) * (-t110 - t159) + g(2) * t125) * t122) * MDP(21) + (-g(1) * (-t94 * pkin(5) - t93 * qJ(6)) - g(2) * (t96 * pkin(5) + t95 * qJ(6) + t103) + (-g(1) * t137 - g(2) * t131) * t124 + (-g(1) * (-t110 - t131) - g(2) * t137) * t122) * MDP(32) + t158 * (g(1) * t94 - g(2) * t96) - t157 * (g(1) * t93 - g(2) * t95) + (-t121 * MDP(10) + t114 * MDP(16) + t123 * MDP(9) - t113 * t160 + MDP(2)) * (g(1) * t122 - g(2) * t124); (-g(3) * (t115 + t138 - t148) + t134 * (t130 * t113 + t114 * t120 + t156)) * MDP(32) + (g(3) * t121 + t134 * t123) * MDP(10) + (-g(3) * t123 + t134 * t121) * MDP(9) + (-g(1) * (t135 * t124 + t102) - g(2) * (t135 * t122 + t101) - g(3) * (t115 + t159)) * MDP(21) + t127; (-g(1) * (-t124 * t155 + t102) - g(2) * (-t122 * t155 + t101) - g(3) * t159) * MDP(21) + (-g(3) * t138 + t120 * t129 + (g(3) * t120 + t134 * t130) * t113) * MDP(32) + t127; (-MDP(21) - MDP(32)) * t91; (-g(1) * (-pkin(5) * t95 + qJ(6) * t96) - g(2) * (-pkin(5) * t93 + qJ(6) * t94) - (-pkin(5) * t111 + qJ(6) * t112) * t151) * MDP(32) + t158 * t77 + t157 * (g(1) * t96 + g(2) * t94 + t112 * t151); -t77 * MDP(32);];
taug  = t1;
