% Calculate Gravitation load on the joints for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:43
% EndTime: 2019-03-09 10:47:45
% DurationCPUTime: 0.43s
% Computational Cost: add. (211->83), mult. (393->126), div. (0->0), fcn. (392->10), ass. (0->51)
t115 = qJ(4) + pkin(10);
t107 = sin(t115);
t108 = cos(t115);
t117 = sin(qJ(6));
t121 = cos(qJ(6));
t118 = sin(qJ(4));
t119 = sin(qJ(2));
t122 = cos(qJ(4));
t123 = cos(qJ(2));
t129 = t118 * t123 - t119 * t122;
t130 = t119 * t107 + t123 * t108;
t124 = cos(qJ(1));
t139 = t123 * t124;
t120 = sin(qJ(1));
t140 = t120 * t123;
t141 = t119 * t124;
t142 = t119 * t120;
t143 = t118 * t119;
t93 = t122 * t123 + t143;
t150 = g(3) * t93;
t85 = t129 * t120;
t86 = t93 * t120;
t134 = t118 * t139;
t87 = -t122 * t141 + t134;
t88 = t93 * t124;
t156 = (MDP(29) * t121 - MDP(30) * t117) * (g(1) * (t107 * t139 - t108 * t141) + g(2) * (t107 * t140 - t108 * t142) + g(3) * t130) + (g(1) * t87 + g(2) * t85 + t150) * MDP(20) + (g(1) * t88 + g(2) * t86 - g(3) * t129) * MDP(21);
t138 = t123 * pkin(2) + t119 * qJ(3);
t98 = g(1) * t124 + g(2) * t120;
t154 = MDP(9) + MDP(11);
t153 = MDP(10) - MDP(13);
t151 = g(3) * (-t107 * t123 + t108 * t119);
t149 = pkin(4) * t118;
t148 = g(1) * t120;
t106 = pkin(4) * t122 + pkin(3);
t144 = t106 * t123;
t135 = pkin(4) * t143;
t133 = t124 * pkin(1) + pkin(2) * t139 + t120 * pkin(7) + qJ(3) * t141;
t97 = -g(2) * t124 + t148;
t80 = t130 * t120;
t132 = t80 * t117 - t121 * t124;
t131 = t117 * t124 + t80 * t121;
t128 = -pkin(1) - t138;
t116 = -qJ(5) - pkin(8);
t112 = t124 * pkin(7);
t104 = qJ(3) * t139;
t102 = qJ(3) * t140;
t83 = -g(3) * t123 + t98 * t119;
t82 = t130 * t124;
t78 = -t117 * t120 + t121 * t82;
t77 = -t117 * t82 - t120 * t121;
t1 = [(-g(1) * t112 - g(2) * t133 - t128 * t148) * MDP(14) + (g(1) * t86 - g(2) * t88) * MDP(20) + (-g(1) * t85 + g(2) * t87) * MDP(21) + (-g(1) * (t116 * t124 + t112) - g(2) * (t106 * t139 + t124 * t135 + t133) + (-g(1) * (t128 - t135 - t144) - g(2) * t116) * t120) * MDP(23) + (g(1) * t131 - g(2) * t78) * MDP(29) + (-g(1) * t132 - g(2) * t77) * MDP(30) + (MDP(3) - MDP(12) + MDP(22)) * t98 + (-t153 * t119 + t154 * t123 + MDP(2)) * t97; (-g(1) * (-pkin(2) * t141 + t104) - g(2) * (-pkin(2) * t142 + t102) - g(3) * t138) * MDP(14) + (-g(1) * (pkin(4) * t134 + t104) - g(2) * (t140 * t149 + t102) - g(3) * (t138 + t144) + (-g(3) * t149 + t98 * (pkin(2) + t106)) * t119) * MDP(23) + t153 * (g(3) * t119 + t98 * t123) + t154 * t83 - t156; (-MDP(14) - MDP(23)) * t83; (t98 * t129 + t150) * pkin(4) * MDP(23) + t156; t97 * MDP(23); (-g(1) * t77 + g(2) * t132 + t117 * t151) * MDP(29) + (g(1) * t78 + g(2) * t131 + t121 * t151) * MDP(30);];
taug  = t1;
