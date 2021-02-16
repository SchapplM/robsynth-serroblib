% Calculate Gravitation load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:06
% EndTime: 2021-01-16 02:06:11
% DurationCPUTime: 0.95s
% Computational Cost: add. (396->110), mult. (668->179), div. (0->0), fcn. (767->14), ass. (0->57)
t115 = sin(qJ(3));
t117 = cos(qJ(3));
t112 = sin(pkin(6));
t116 = sin(qJ(2));
t139 = t112 * t116;
t147 = cos(pkin(6));
t157 = -t115 * t139 + t147 * t117;
t111 = sin(pkin(10));
t138 = t112 * t117;
t131 = t111 * t147;
t146 = cos(pkin(10));
t152 = cos(qJ(2));
t90 = t116 * t131 - t146 * t152;
t156 = t111 * t138 + t90 * t115;
t109 = qJ(3) + pkin(11);
t105 = sin(t109);
t107 = cos(t109);
t130 = t112 * t146;
t126 = t147 * t146;
t92 = t111 * t152 + t116 * t126;
t81 = t92 * t105 + t107 * t130;
t140 = t111 * t112;
t83 = t105 * t90 + t107 * t140;
t87 = t105 * t139 - t147 * t107;
t155 = g(1) * t83 - g(2) * t81 - g(3) * t87;
t154 = -MDP(13) + MDP(18);
t151 = g(3) * t112;
t103 = pkin(3) * t117 + pkin(2);
t114 = qJ(4) + pkin(8);
t93 = t146 * t116 + t152 * t131;
t150 = -t93 * t103 - t90 * t114;
t133 = t112 * t152;
t149 = t103 * t133 + t114 * t139;
t145 = qJ(5) * t105;
t108 = pkin(12) + qJ(6);
t104 = sin(t108);
t144 = t104 * t107;
t106 = cos(t108);
t143 = t106 * t107;
t110 = sin(pkin(12));
t142 = t107 * t110;
t113 = cos(pkin(12));
t141 = t107 * t113;
t137 = MDP(15) + MDP(19);
t134 = t107 * t152;
t80 = -t105 * t140 + t107 * t90;
t91 = t111 * t116 - t152 * t126;
t132 = -t91 * t103 + t114 * t92;
t128 = t156 * pkin(3);
t125 = -pkin(4) * t107 - t145;
t124 = t157 * pkin(3);
t121 = -t92 * t115 - t117 * t130;
t119 = t121 * pkin(3);
t118 = g(1) * t93 + g(2) * t91 - g(3) * t133;
t88 = t147 * t105 + t107 * t139;
t82 = -t105 * t130 + t107 * t92;
t1 = [(-MDP(1) - t137) * g(3); (-g(1) * t150 - g(2) * t132 - g(3) * t149) * MDP(15) + (-g(1) * (-t110 * t90 - t93 * t141) - g(2) * (t110 * t92 - t91 * t141) - (t110 * t116 + t113 * t134) * t151) * MDP(16) + (-g(1) * (-t113 * t90 + t93 * t142) - g(2) * (t113 * t92 + t91 * t142) - (-t110 * t134 + t113 * t116) * t151) * MDP(17) + (-g(1) * (t125 * t93 + t150) - g(2) * (t125 * t91 + t132) - g(3) * ((pkin(4) * t134 + t152 * t145) * t112 + t149)) * MDP(19) + (-g(1) * (-t104 * t90 - t93 * t143) - g(2) * (t104 * t92 - t91 * t143) - (t104 * t116 + t106 * t134) * t151) * MDP(25) + (-g(1) * (-t106 * t90 + t93 * t144) - g(2) * (t106 * t92 + t91 * t144) - (-t104 * t134 + t106 * t116) * t151) * MDP(26) + (MDP(4) - MDP(14)) * (-g(1) * t90 + g(2) * t92 + g(3) * t139) + (t117 * MDP(10) - t115 * MDP(11) + t107 * MDP(12) + t154 * t105 + MDP(3)) * t118; (-g(1) * t156 - g(2) * t121 - g(3) * t157) * MDP(10) + (-g(1) * (-t115 * t140 + t117 * t90) - g(2) * (t115 * t130 - t92 * t117) - g(3) * (-t147 * t115 - t116 * t138)) * MDP(11) + (-g(1) * t128 - g(2) * t119 - g(3) * t124) * MDP(15) + (-g(1) * (pkin(4) * t83 - qJ(5) * t80 + t128) - g(2) * (-t81 * pkin(4) + t82 * qJ(5) + t119) - g(3) * (-pkin(4) * t87 + qJ(5) * t88 + t124)) * MDP(19) + t154 * (g(1) * t80 - g(2) * t82 - g(3) * t88) + (-MDP(16) * t113 + MDP(17) * t110 - MDP(25) * t106 + MDP(26) * t104 - MDP(12)) * t155; -t137 * t118; t155 * MDP(19); (-g(1) * (t104 * t80 + t106 * t93) - g(2) * (-t104 * t82 + t106 * t91) - g(3) * (-t88 * t104 - t106 * t133)) * MDP(25) + (-g(1) * (-t104 * t93 + t106 * t80) - g(2) * (-t104 * t91 - t106 * t82) - g(3) * (t104 * t133 - t88 * t106)) * MDP(26);];
taug = t1;
