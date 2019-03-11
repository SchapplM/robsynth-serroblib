% Calculate Gravitation load on the joints for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:07
% EndTime: 2019-03-08 20:03:09
% DurationCPUTime: 0.64s
% Computational Cost: add. (517->92), mult. (1342->153), div. (0->0), fcn. (1705->12), ass. (0->56)
t181 = MDP(12) - MDP(21);
t179 = MDP(18) + MDP(20);
t178 = MDP(19) - MDP(22);
t135 = sin(pkin(11));
t142 = sin(qJ(2));
t145 = cos(qJ(2));
t170 = cos(pkin(11));
t153 = -t142 * t135 + t145 * t170;
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t139 = cos(pkin(6));
t163 = t139 * t145;
t180 = -t136 * t163 - t138 * t142;
t127 = t153 * t139;
t130 = -t145 * t135 - t142 * t170;
t115 = t127 * t138 + t130 * t136;
t118 = -t127 * t136 + t130 * t138;
t137 = sin(pkin(6));
t125 = t153 * t137;
t177 = -g(1) * t118 - g(2) * t115 - g(3) * t125;
t169 = t136 * t142;
t141 = sin(qJ(4));
t168 = t137 * t141;
t144 = cos(qJ(4));
t167 = t137 * t144;
t166 = t137 * t145;
t164 = t139 * t142;
t140 = sin(qJ(5));
t162 = t140 * t144;
t143 = cos(qJ(5));
t160 = t143 * t144;
t159 = MDP(23) + MDP(5);
t157 = t138 * t163;
t128 = t130 * t139;
t116 = -t128 * t138 + t136 * t153;
t117 = -t128 * t136 - t138 * t153;
t126 = t130 * t137;
t121 = -t126 * t144 + t139 * t141;
t102 = t121 * t140 + t125 * t143;
t105 = t116 * t144 - t138 * t168;
t94 = t105 * t140 + t115 * t143;
t107 = -t117 * t144 + t136 * t168;
t96 = t107 * t140 + t118 * t143;
t152 = g(1) * t96 + g(2) * t94 + g(3) * t102;
t146 = -g(1) * t180 - g(3) * t166;
t131 = pkin(2) * t157;
t109 = t125 * t160 - t126 * t140;
t108 = t125 * t162 + t126 * t143;
t103 = t121 * t143 - t125 * t140;
t101 = -t117 * t140 + t118 * t160;
t100 = t117 * t143 + t118 * t162;
t99 = t115 * t160 + t116 * t140;
t98 = t115 * t162 - t116 * t143;
t97 = t107 * t143 - t118 * t140;
t95 = t105 * t143 - t115 * t140;
t1 = [(-MDP(1) - t159) * g(3); (-g(2) * (t157 - t169) + t146) * MDP(3) + (-g(1) * (t136 * t164 - t138 * t145) - g(2) * (-t136 * t145 - t138 * t164) + g(3) * t137 * t142) * MDP(4) + (-g(2) * t131 + (g(2) * t169 + t146) * pkin(2)) * MDP(5) + (-g(1) * (pkin(2) * t180 + t101 * pkin(5) - t117 * pkin(8) + t100 * qJ(6)) - g(2) * (-pkin(2) * t169 + t99 * pkin(5) + pkin(8) * t116 + t98 * qJ(6) + t131) - g(3) * (pkin(2) * t166 + t109 * pkin(5) - t126 * pkin(8) + t108 * qJ(6)) + t177 * (pkin(4) * t144 + pkin(9) * t141 + pkin(3))) * MDP(23) + t179 * (-g(1) * t101 - g(2) * t99 - g(3) * t109) + t178 * (g(1) * t100 + g(2) * t98 + g(3) * t108) + (t144 * MDP(11) - t181 * t141) * t177; t159 * (-g(3) * t139 + (-g(1) * t136 + g(2) * t138) * t137); (-pkin(9) * MDP(23) + t181) * (g(1) * t107 + g(2) * t105 + g(3) * t121) + (-MDP(11) - t179 * t143 + t178 * t140 - MDP(23) * (pkin(5) * t143 + qJ(6) * t140 + pkin(4))) * (g(3) * (t126 * t141 + t139 * t144) + g(2) * (-t116 * t141 - t138 * t167) + g(1) * (t117 * t141 + t136 * t167)); (-g(1) * (-pkin(5) * t96 + qJ(6) * t97) - g(2) * (-pkin(5) * t94 + qJ(6) * t95) - g(3) * (-t102 * pkin(5) + qJ(6) * t103)) * MDP(23) + t179 * t152 + t178 * (g(1) * t97 + g(2) * t95 + g(3) * t103); -t152 * MDP(23);];
taug  = t1;
