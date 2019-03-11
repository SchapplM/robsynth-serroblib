% Calculate Gravitation load on the joints for
% S6PRPRRP4
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:19
% EndTime: 2019-03-08 20:12:21
% DurationCPUTime: 0.62s
% Computational Cost: add. (463->88), mult. (818->139), div. (0->0), fcn. (972->12), ass. (0->51)
t158 = MDP(15) - MDP(24);
t156 = MDP(21) + MDP(23);
t155 = MDP(22) - MDP(25);
t119 = sin(pkin(10));
t124 = sin(qJ(2));
t126 = cos(qJ(2));
t147 = cos(pkin(10));
t148 = cos(pkin(6));
t136 = t148 * t147;
t105 = t119 * t124 - t126 * t136;
t138 = t119 * t148;
t107 = t147 * t124 + t126 * t138;
t157 = -g(1) * t107 - g(2) * t105;
t120 = sin(pkin(6));
t149 = g(3) * t120;
t117 = pkin(11) + qJ(4);
t116 = cos(t117);
t123 = sin(qJ(5));
t146 = t116 * t123;
t125 = cos(qJ(5));
t145 = t116 * t125;
t144 = t119 * t120;
t143 = t120 * t124;
t142 = t120 * t126;
t141 = t125 * t126;
t140 = MDP(26) + MDP(8);
t139 = t123 * t142;
t137 = t120 * t147;
t106 = t119 * t126 + t124 * t136;
t115 = sin(t117);
t92 = t106 * t116 - t115 * t137;
t82 = -t105 * t125 + t123 * t92;
t108 = -t124 * t138 + t147 * t126;
t94 = t108 * t116 + t115 * t144;
t84 = -t107 * t125 + t123 * t94;
t100 = t148 * t115 + t116 * t143;
t95 = t100 * t123 + t120 * t141;
t133 = g(1) * t84 + g(2) * t82 + g(3) * t95;
t128 = g(3) * t142 + t157;
t122 = -pkin(8) - qJ(3);
t121 = cos(pkin(11));
t98 = (t116 * t141 + t123 * t124) * t120;
t97 = t116 * t139 - t125 * t143;
t96 = t100 * t125 - t139;
t89 = -t107 * t145 + t108 * t123;
t88 = -t107 * t146 - t108 * t125;
t87 = -t105 * t145 + t106 * t123;
t86 = -t105 * t146 - t106 * t125;
t85 = t107 * t123 + t125 * t94;
t83 = t105 * t123 + t125 * t92;
t1 = [(-MDP(1) - t140) * g(3); (-g(1) * (-pkin(2) * t107 + qJ(3) * t108) - g(2) * (-pkin(2) * t105 + qJ(3) * t106) - (pkin(2) * t126 + qJ(3) * t124) * t149) * MDP(8) + (-g(1) * (pkin(5) * t89 + qJ(6) * t88 - t108 * t122) - g(2) * (pkin(5) * t87 + qJ(6) * t86 - t106 * t122) - g(3) * (t98 * pkin(5) + t97 * qJ(6)) + t122 * t124 * t149 + (-t126 * t149 - t157) * (pkin(3) * t121 + pkin(4) * t116 + pkin(9) * t115 + pkin(2))) * MDP(26) + t156 * (-g(1) * t89 - g(2) * t87 - g(3) * t98) + t155 * (g(1) * t88 + g(2) * t86 + g(3) * t97) + (MDP(4) - MDP(7)) * (g(1) * t108 + g(2) * t106 + g(3) * t143) + (-t116 * MDP(14) - t121 * MDP(5) + MDP(6) * sin(pkin(11)) + t158 * t115 - MDP(3)) * t128; t140 * t128; (-pkin(9) * MDP(26) + t158) * (g(1) * t94 + g(2) * t92 + g(3) * t100) + (-MDP(14) - t156 * t125 + t155 * t123 - MDP(26) * (pkin(5) * t125 + qJ(6) * t123 + pkin(4))) * (g(3) * (-t115 * t143 + t148 * t116) + g(2) * (-t106 * t115 - t116 * t137) + g(1) * (-t108 * t115 + t116 * t144)); (-g(1) * (-pkin(5) * t84 + qJ(6) * t85) - g(2) * (-pkin(5) * t82 + qJ(6) * t83) - g(3) * (-pkin(5) * t95 + qJ(6) * t96)) * MDP(26) + t156 * t133 + t155 * (g(1) * t85 + g(2) * t83 + g(3) * t96); -t133 * MDP(26);];
taug  = t1;
