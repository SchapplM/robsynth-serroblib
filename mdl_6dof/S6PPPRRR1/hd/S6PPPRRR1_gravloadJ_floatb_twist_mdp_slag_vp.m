% Calculate Gravitation load on the joints for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:35
% EndTime: 2019-03-08 18:40:37
% DurationCPUTime: 0.36s
% Computational Cost: add. (681->78), mult. (1944->147), div. (0->0), fcn. (2565->18), ass. (0->60)
t106 = sin(pkin(14));
t147 = sin(pkin(12));
t151 = cos(pkin(13));
t136 = t147 * t151;
t146 = sin(pkin(13));
t152 = cos(pkin(12));
t137 = t152 * t146;
t155 = cos(pkin(6));
t127 = t155 * t137 + t136;
t150 = cos(pkin(14));
t135 = t147 * t146;
t138 = t152 * t151;
t126 = t155 * t138 - t135;
t107 = sin(pkin(6));
t149 = sin(pkin(7));
t142 = t107 * t149;
t154 = cos(pkin(7));
t158 = t126 * t154 - t152 * t142;
t115 = t127 * t106 - t158 * t150;
t120 = t152 * t107 * t154 + t126 * t149;
t148 = sin(pkin(8));
t153 = cos(pkin(8));
t161 = t115 * t153 + t120 * t148;
t129 = -t155 * t135 + t138;
t128 = -t155 * t136 - t137;
t141 = t107 * t147;
t157 = t128 * t154 + t149 * t141;
t116 = t129 * t106 - t157 * t150;
t121 = t128 * t149 - t154 * t141;
t160 = t116 * t153 + t121 * t148;
t139 = t154 * t151;
t140 = t155 * t149;
t122 = t150 * t140 + (-t146 * t106 + t150 * t139) * t107;
t130 = t151 * t142 - t155 * t154;
t159 = t122 * t153 - t130 * t148;
t156 = cos(qJ(4));
t108 = sin(qJ(6));
t112 = cos(qJ(5));
t145 = t108 * t112;
t111 = cos(qJ(6));
t144 = t111 * t112;
t143 = MDP(2) + MDP(3);
t110 = sin(qJ(4));
t109 = sin(qJ(5));
t103 = t107 * t146 * t150 + (t107 * t139 + t140) * t106;
t99 = -t122 * t148 - t130 * t153;
t98 = t157 * t106 + t129 * t150;
t97 = t158 * t106 + t127 * t150;
t94 = t116 * t148 - t121 * t153;
t93 = t115 * t148 - t120 * t153;
t92 = t103 * t156 + t159 * t110;
t91 = t103 * t110 - t159 * t156;
t90 = -t160 * t110 + t98 * t156;
t89 = t98 * t110 + t160 * t156;
t88 = -t161 * t110 + t97 * t156;
t87 = t97 * t110 + t161 * t156;
t86 = t99 * t109 + t92 * t112;
t84 = t94 * t109 + t90 * t112;
t82 = t93 * t109 + t88 * t112;
t1 = [(-MDP(1) - t143) * g(3); t143 * (-g(3) * t155 + (-t147 * g(1) + t152 * g(2)) * t107); (g(1) * t121 + g(2) * t120 + g(3) * t130) * MDP(3); (g(1) * t90 + g(2) * t88 + g(3) * t92) * MDP(6) + (-g(1) * (t90 * t108 - t89 * t144) - g(2) * (t88 * t108 - t87 * t144) - g(3) * (t92 * t108 - t91 * t144)) * MDP(19) + (-g(1) * (t90 * t111 + t89 * t145) - g(2) * (t88 * t111 + t87 * t145) - g(3) * (t92 * t111 + t91 * t145)) * MDP(20) + (t112 * MDP(12) - MDP(13) * t109 + MDP(5)) * (g(1) * t89 + g(2) * t87 + g(3) * t91); (g(1) * t84 + g(2) * t82 + g(3) * t86) * MDP(13) + (-MDP(19) * t111 + MDP(20) * t108 - MDP(12)) * (g(1) * (-t90 * t109 + t94 * t112) + g(2) * (-t88 * t109 + t93 * t112) + g(3) * (-t92 * t109 + t99 * t112)); (-g(1) * (-t84 * t108 + t89 * t111) - g(2) * (-t82 * t108 + t87 * t111) - g(3) * (-t86 * t108 + t91 * t111)) * MDP(19) + (-g(1) * (-t89 * t108 - t84 * t111) - g(2) * (-t87 * t108 - t82 * t111) - g(3) * (-t91 * t108 - t86 * t111)) * MDP(20);];
taug  = t1;
