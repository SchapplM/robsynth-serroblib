% Calculate Gravitation load on the joints for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:47
% EndTime: 2019-03-08 23:36:49
% DurationCPUTime: 0.79s
% Computational Cost: add. (456->103), mult. (1197->171), div. (0->0), fcn. (1491->12), ass. (0->49)
t167 = MDP(11) - MDP(20);
t164 = MDP(17) + MDP(19);
t163 = MDP(18) - MDP(21);
t122 = sin(qJ(3));
t123 = sin(qJ(2));
t119 = sin(pkin(6));
t126 = cos(qJ(3));
t152 = t119 * t126;
t155 = cos(pkin(6));
t112 = t155 * t122 + t123 * t152;
t125 = cos(qJ(4));
t121 = sin(qJ(4));
t127 = cos(qJ(2));
t151 = t119 * t127;
t147 = t121 * t151;
t100 = t112 * t125 - t147;
t120 = sin(qJ(6));
t124 = cos(qJ(6));
t118 = sin(pkin(11));
t154 = cos(pkin(11));
t144 = t155 * t154;
t107 = t118 * t123 - t127 * t144;
t108 = t118 * t127 + t123 * t144;
t145 = t119 * t154;
t96 = t108 * t126 - t122 * t145;
t86 = -t107 * t125 + t121 * t96;
t87 = t107 * t121 + t125 * t96;
t146 = t118 * t155;
t109 = t154 * t123 + t127 * t146;
t110 = -t123 * t146 + t154 * t127;
t98 = t118 * t119 * t122 + t110 * t126;
t88 = -t109 * t125 + t121 * t98;
t89 = t109 * t121 + t125 * t98;
t148 = t125 * t127;
t99 = t112 * t121 + t119 * t148;
t166 = (g(1) * (t120 * t89 - t124 * t88) + g(2) * (t120 * t87 - t124 * t86) + g(3) * (t100 * t120 - t124 * t99)) * MDP(28) + (g(1) * (t120 * t88 + t124 * t89) + g(2) * (t120 * t86 + t124 * t87) + g(3) * (t100 * t124 + t120 * t99)) * MDP(29);
t165 = g(1) * t109 + g(2) * t107;
t153 = t119 * t123;
t150 = t121 * t126;
t149 = t125 * t126;
t135 = pkin(3) * t126 + pkin(9) * t122 + pkin(2);
t133 = g(1) * t88 + g(2) * t86 + g(3) * t99;
t102 = (t121 * t123 + t126 * t148) * t119;
t101 = -t125 * t153 + t126 * t147;
t94 = -t109 * t149 + t110 * t121;
t93 = -t109 * t150 - t110 * t125;
t92 = -t107 * t149 + t108 * t121;
t91 = -t107 * t150 - t108 * t125;
t1 = [(-MDP(1) - MDP(22)) * g(3); (g(1) * t110 + g(2) * t108 + g(3) * t153) * MDP(4) + (-g(1) * (t94 * pkin(4) + t110 * pkin(8) + t93 * qJ(5)) - g(2) * (t92 * pkin(4) + pkin(8) * t108 + qJ(5) * t91) + t165 * t135 + (-t102 * pkin(4) - t101 * qJ(5) - (pkin(8) * t123 + t135 * t127) * t119) * g(3)) * MDP(22) + (-g(1) * (t120 * t93 + t124 * t94) - g(2) * (t120 * t91 + t124 * t92) - g(3) * (t101 * t120 + t102 * t124)) * MDP(28) + (-g(1) * (-t120 * t94 + t124 * t93) - g(2) * (-t120 * t92 + t124 * t91) - g(3) * (t101 * t124 - t102 * t120)) * MDP(29) + t164 * (-g(1) * t94 - g(2) * t92 - g(3) * t102) + t163 * (g(1) * t93 + g(2) * t91 + g(3) * t101) + (-t126 * MDP(10) + t122 * t167 - MDP(3)) * (g(3) * t151 - t165); (-pkin(9) * MDP(22) + t167) * (g(1) * t98 + g(2) * t96 + g(3) * t112) + (-MDP(10) - t164 * t125 + t163 * t121 - (pkin(4) * t125 + qJ(5) * t121 + pkin(3)) * MDP(22) + MDP(29) * (t120 * t125 - t121 * t124) - (t120 * t121 + t124 * t125) * MDP(28)) * (g(3) * (-t122 * t153 + t155 * t126) + g(2) * (-t108 * t122 - t126 * t145) + g(1) * (-t110 * t122 + t118 * t152)); (-g(1) * (-pkin(4) * t88 + qJ(5) * t89) - g(2) * (-pkin(4) * t86 + qJ(5) * t87) - g(3) * (-pkin(4) * t99 + qJ(5) * t100)) * MDP(22) + t164 * t133 + t163 * (g(1) * t89 + g(2) * t87 + g(3) * t100) - t166; -t133 * MDP(22); t166;];
taug  = t1;
