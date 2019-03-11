% Calculate Gravitation load on the joints for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:18
% EndTime: 2019-03-09 10:06:20
% DurationCPUTime: 0.68s
% Computational Cost: add. (245->98), mult. (551->132), div. (0->0), fcn. (531->6), ass. (0->52)
t162 = MDP(21) - MDP(24) - MDP(27);
t120 = sin(qJ(1));
t123 = cos(qJ(1));
t100 = g(1) * t123 + g(2) * t120;
t161 = MDP(10) - MDP(13);
t159 = MDP(20) + MDP(22) + MDP(26);
t119 = sin(qJ(2));
t122 = cos(qJ(2));
t87 = g(3) * t119 + t100 * t122;
t157 = MDP(9) - MDP(12) + MDP(23) - MDP(28);
t156 = pkin(2) + pkin(8);
t155 = g(1) * t120;
t151 = g(3) * t122;
t113 = t122 * pkin(2);
t142 = t120 * t122;
t102 = qJ(3) * t142;
t118 = sin(qJ(4));
t146 = t118 * t122;
t135 = pkin(4) * t146;
t150 = t120 * t135 + t102;
t139 = t122 * t123;
t104 = qJ(3) * t139;
t149 = t123 * t135 + t104;
t148 = qJ(5) * t118;
t121 = cos(qJ(4));
t147 = qJ(5) * t121;
t109 = t119 * qJ(3);
t145 = t119 * t120;
t144 = t119 * t123;
t143 = t120 * t121;
t141 = t121 * t122;
t140 = t121 * t123;
t138 = t113 + t109;
t137 = MDP(25) + MDP(29);
t136 = -qJ(6) + t156;
t134 = qJ(5) * t141;
t133 = -pkin(1) - t109;
t92 = t118 * t120 - t119 * t140;
t93 = t118 * t144 + t143;
t132 = -t92 * pkin(4) + qJ(5) * t93;
t94 = t118 * t123 + t119 * t143;
t95 = -t118 * t145 + t140;
t131 = t94 * pkin(4) - qJ(5) * t95;
t130 = t123 * pkin(1) + pkin(2) * t139 + t120 * pkin(7) + qJ(3) * t144;
t127 = pkin(5) * t118 - t147;
t126 = g(3) * (t119 * t118 * pkin(4) + t122 * pkin(8) + t138);
t114 = t123 * pkin(7);
t125 = t123 * pkin(3) + t95 * pkin(4) + qJ(5) * t94 + t114;
t80 = g(1) * t92 - g(2) * t94 + g(3) * t141;
t124 = t120 * pkin(3) + t93 * pkin(4) + pkin(8) * t139 + t92 * qJ(5) + t130;
t86 = t100 * t119 - t151;
t1 = [(-g(1) * t114 - g(2) * t130 - (t133 - t113) * t155) * MDP(14) + (-g(1) * t125 - g(2) * t124 - (-t156 * t122 + t133) * t155) * MDP(25) + (-g(1) * (pkin(5) * t95 + t125) - g(2) * (t93 * pkin(5) - qJ(6) * t139 + t124) - (-t136 * t122 + t133) * t155) * MDP(29) + (MDP(3) - MDP(11)) * t100 + t159 * (-g(1) * t95 - g(2) * t93) + t157 * (g(1) * t142 - g(2) * t139) + t162 * (g(1) * t94 + g(2) * t92) + (-t119 * t161 + MDP(2)) * (-g(2) * t123 + t155); (-g(1) * (-pkin(2) * t144 + t104) - g(2) * (-pkin(2) * t145 + t102) - g(3) * t138) * MDP(14) + (-g(1) * (-t123 * t134 + t149) - g(2) * (-t120 * t134 + t150) - t126 + (g(3) * t147 + t100 * t156) * t119) * MDP(25) + (-g(1) * t149 - g(2) * t150 - t126 + (g(3) * qJ(6) - t100 * t127) * t122 + (-g(3) * t127 + t100 * t136) * t119) * MDP(29) + t157 * t86 + (-t118 * t159 - t121 * t162 + t161) * t87; (-MDP(14) - t137) * t86; (-g(1) * t132 - g(2) * t131 - (-pkin(4) * t121 - t148) * t151) * MDP(25) + (-g(1) * (-pkin(5) * t92 + t132) - g(2) * (pkin(5) * t94 + t131) - (-t148 + (-pkin(4) - pkin(5)) * t121) * t151) * MDP(29) - t162 * (-g(1) * t93 + g(2) * t95 + g(3) * t146) + t159 * t80; -t137 * t80; t87 * MDP(29);];
taug  = t1;
