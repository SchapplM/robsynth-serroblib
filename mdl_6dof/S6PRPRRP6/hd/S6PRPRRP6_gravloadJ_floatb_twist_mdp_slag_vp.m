% Calculate Gravitation load on the joints for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:56
% EndTime: 2019-03-08 20:20:58
% DurationCPUTime: 0.51s
% Computational Cost: add. (316->87), mult. (808->135), div. (0->0), fcn. (963->10), ass. (0->51)
t159 = -MDP(14) + MDP(23);
t157 = MDP(20) + MDP(22);
t156 = MDP(21) - MDP(24);
t120 = sin(pkin(10));
t122 = cos(pkin(10));
t129 = cos(qJ(2));
t123 = cos(pkin(6));
t126 = sin(qJ(2));
t146 = t123 * t126;
t108 = t120 * t129 + t122 * t146;
t110 = -t120 * t146 + t122 * t129;
t158 = -g(1) * t110 - g(2) * t108;
t121 = sin(pkin(6));
t125 = sin(qJ(4));
t150 = t121 * t125;
t149 = t121 * t126;
t128 = cos(qJ(4));
t148 = t121 * t128;
t147 = t121 * t129;
t145 = t123 * t129;
t124 = sin(qJ(5));
t144 = t124 * t125;
t127 = cos(qJ(5));
t143 = t125 * t127;
t142 = t126 * t127;
t141 = pkin(2) * t147 + qJ(3) * t149;
t140 = MDP(25) + MDP(7);
t139 = t124 * t149;
t138 = pkin(4) * t125 - pkin(9) * t128;
t109 = t120 * t145 + t122 * t126;
t93 = t109 * t125 + t120 * t148;
t82 = -t110 * t127 + t124 * t93;
t107 = t120 * t126 - t122 * t145;
t95 = -t107 * t125 + t122 * t148;
t84 = -t108 * t127 - t124 * t95;
t112 = t123 * t128 - t125 * t147;
t96 = t112 * t124 - t121 * t142;
t135 = g(1) * t82 + g(2) * t84 + g(3) * t96;
t87 = -g(1) * t109 - g(2) * t107 + g(3) * t147;
t106 = t109 * pkin(2);
t105 = t107 * pkin(2);
t99 = (t124 * t129 + t125 * t142) * t121;
t98 = t125 * t139 - t127 * t147;
t97 = t112 * t127 + t139;
t91 = -t109 * t124 + t110 * t143;
t90 = t109 * t127 + t110 * t144;
t89 = -t107 * t124 + t108 * t143;
t88 = t107 * t127 + t108 * t144;
t85 = t108 * t124 - t127 * t95;
t83 = t110 * t124 + t127 * t93;
t1 = [(-MDP(1) - t140) * g(3); (-g(1) * (qJ(3) * t110 - t106) - g(2) * (qJ(3) * t108 - t105) - g(3) * t141) * MDP(7) + (-g(1) * (pkin(5) * t91 - t109 * pkin(8) + t90 * qJ(6) - t106) - g(2) * (pkin(5) * t89 - pkin(8) * t107 + qJ(6) * t88 - t105) + t158 * (qJ(3) + t138) + (-t99 * pkin(5) - t98 * qJ(6) - t141 - (pkin(8) * t129 + t138 * t126) * t121) * g(3)) * MDP(25) + (-MDP(3) + MDP(5)) * t87 + t157 * (-g(1) * t91 - g(2) * t89 - g(3) * t99) + t156 * (g(1) * t90 + g(2) * t88 + g(3) * t98) + (-MDP(13) * t125 + t159 * t128 + MDP(4) - MDP(6)) * (g(3) * t149 - t158); t140 * t87; (-pkin(9) * MDP(25) - t159) * (g(1) * t93 - g(2) * t95 + g(3) * t112) + (-MDP(13) - t157 * t127 + t156 * t124 - MDP(25) * (pkin(5) * t127 + qJ(6) * t124 + pkin(4))) * (g(3) * (-t123 * t125 - t128 * t147) + g(2) * (t107 * t128 + t122 * t150) + g(1) * (t109 * t128 - t120 * t150)); (-g(1) * (-pkin(5) * t82 + qJ(6) * t83) - g(2) * (-pkin(5) * t84 + qJ(6) * t85) - g(3) * (-pkin(5) * t96 + qJ(6) * t97)) * MDP(25) + t157 * t135 + t156 * (g(1) * t83 + g(2) * t85 + g(3) * t97); -t135 * MDP(25);];
taug  = t1;
