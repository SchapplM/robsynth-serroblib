% Calculate Gravitation load on the joints for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:26
% EndTime: 2021-01-15 16:45:32
% DurationCPUTime: 0.71s
% Computational Cost: add. (244->107), mult. (617->172), div. (0->0), fcn. (723->10), ass. (0->61)
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t130 = cos(qJ(2));
t123 = cos(pkin(5));
t127 = sin(qJ(2));
t144 = t123 * t127;
t108 = t120 * t130 + t122 * t144;
t126 = sin(qJ(3));
t121 = sin(pkin(5));
t129 = cos(qJ(3));
t147 = t121 * t129;
t85 = t108 * t126 + t122 * t147;
t106 = t120 * t144 - t122 * t130;
t87 = t106 * t126 + t120 * t147;
t165 = g(1) * t87 - g(2) * t85;
t154 = g(3) * t121;
t141 = t126 * t127;
t111 = t121 * t141 - t123 * t129;
t155 = g(3) * t111;
t164 = -t155 + t165;
t163 = -MDP(11) + MDP(21);
t139 = t127 * t129;
t112 = t121 * t139 + t123 * t126;
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t146 = t121 * t130;
t159 = g(3) * (-t112 * t125 - t128 * t146);
t158 = g(3) * (-t112 * t128 + t125 * t146);
t137 = t129 * t130;
t157 = (-t125 * t137 + t127 * t128) * t154;
t156 = (t125 * t127 + t128 * t137) * t154;
t104 = t106 * t125;
t143 = t123 * t130;
t107 = t120 * t127 - t122 * t143;
t151 = t107 * t125;
t100 = t107 * t128;
t101 = t108 * t125;
t109 = t120 * t143 + t122 * t127;
t150 = t109 * t125;
t103 = t109 * t128;
t149 = t120 * t123;
t148 = t121 * t126;
t145 = t122 * t123;
t142 = t125 * t129;
t140 = t126 * t130;
t138 = t128 * t129;
t84 = t106 * t129 - t120 * t148;
t110 = -t123 * t141 - t147;
t136 = -g(1) * (t120 * t110 + t122 * t140) - g(2) * (-t122 * t110 + t120 * t140);
t119 = pkin(4) * t128 + pkin(3);
t124 = -qJ(5) - pkin(8);
t134 = -t129 * t119 + t126 * t124;
t113 = t123 * t139 - t148;
t105 = t106 * t128;
t102 = t108 * t128;
t96 = t109 * t129;
t95 = t107 * t129;
t92 = -t120 * t113 + t122 * t137;
t90 = t122 * t113 + t120 * t137;
t86 = t108 * t129 - t122 * t148;
t1 = [(-MDP(1) - MDP(22)) * g(3); (-g(1) * t106 + g(2) * t108 + t127 * t154) * MDP(4) + (-g(1) * (-t128 * t96 - t104) - g(2) * (-t128 * t95 + t101) - t156) * MDP(17) + (-g(1) * (t125 * t96 - t105) - g(2) * (t125 * t95 + t102) - t157) * MDP(18) + (-g(1) * (-t109 * t138 - t104) - g(2) * (-t107 * t138 + t101) - t156) * MDP(19) + (-g(1) * (t109 * t142 - t105) - g(2) * (t107 * t142 + t102) - t157) * MDP(20) + (-g(1) * (-pkin(4) * t104 - (pkin(2) * t122 + pkin(7) * t149) * t127 + (-pkin(2) * t149 + pkin(7) * t122) * t130 + t134 * t109) - g(2) * (pkin(4) * t101 - (pkin(2) * t120 - pkin(7) * t145) * t127 + (pkin(2) * t145 + pkin(7) * t120) * t130 + t134 * t107) - ((t125 * pkin(4) + pkin(7)) * t127 + (pkin(2) - t134) * t130) * t154) * MDP(22) + (t129 * MDP(10) + t163 * t126 + MDP(3)) * (g(1) * t109 + g(2) * t107 - g(3) * t146); -t164 * MDP(10) + (-g(1) * (t119 * t87 + t124 * t84) - g(2) * (-t119 * t85 - t124 * t86) - g(3) * (-t111 * t119 - t112 * t124)) * MDP(22) + ((-t136 + t155) * MDP(17) - t164 * MDP(19)) * t128 + t163 * (g(1) * t84 - g(2) * t86 - g(3) * t112) + (-(MDP(18) + MDP(20)) * t155 + t136 * MDP(18) + t165 * MDP(20)) * t125; (-g(1) * (-t125 * t92 + t103) - g(2) * (-t125 * t90 + t100) - t159) * MDP(17) + (-g(1) * (-t128 * t92 - t150) - g(2) * (-t128 * t90 - t151) - t158) * MDP(18) + (-g(1) * (t128 * t84 - t150) - g(2) * (-t128 * t86 - t151) - t158) * MDP(20) + (pkin(4) * MDP(22) + MDP(19)) * (-g(1) * (t125 * t84 + t103) - g(2) * (-t125 * t86 + t100) - t159); t164 * MDP(22);];
taug = t1;
