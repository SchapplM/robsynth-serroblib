% Calculate Gravitation load on the joints for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:18
% EndTime: 2019-03-09 22:16:20
% DurationCPUTime: 0.53s
% Computational Cost: add. (399->83), mult. (578->119), div. (0->0), fcn. (599->10), ass. (0->46)
t118 = sin(qJ(6));
t119 = sin(qJ(4));
t122 = cos(qJ(6));
t123 = cos(qJ(4));
t133 = t118 * t119 + t122 * t123;
t134 = t118 * t123 - t119 * t122;
t117 = qJ(2) + qJ(3);
t115 = cos(t117);
t125 = cos(qJ(1));
t142 = t123 * t125;
t121 = sin(qJ(1));
t144 = t119 * t121;
t101 = t115 * t144 + t142;
t141 = t125 * t119;
t143 = t121 * t123;
t102 = t115 * t143 - t141;
t135 = t101 * t118 + t102 * t122;
t140 = -t101 * t122 + t102 * t118;
t114 = sin(t117);
t148 = g(3) * t114;
t103 = t115 * t141 - t143;
t104 = t115 * t142 + t144;
t84 = t103 * t122 - t104 * t118;
t85 = t103 * t118 + t104 * t122;
t170 = (-g(1) * t84 + g(2) * t140 + t134 * t148) * MDP(34) + (g(1) * t85 + g(2) * t135 + t133 * t148) * MDP(35);
t169 = pkin(4) * t123 + qJ(5) * t119;
t165 = MDP(17) - MDP(26);
t161 = MDP(23) + MDP(25);
t160 = MDP(24) - MDP(27);
t137 = g(1) * t125 + g(2) * t121;
t159 = t137 * t114;
t163 = t115 * pkin(3) + t114 * pkin(9);
t120 = sin(qJ(2));
t153 = pkin(2) * t120;
t151 = pkin(9) * t115;
t139 = t169 * t115 + t163;
t124 = cos(qJ(2));
t116 = t124 * pkin(2);
t131 = t116 + pkin(1) + t163;
t130 = t165 * (t137 * t115 + t148) + (-t133 * MDP(34) + t134 * MDP(35) + t160 * t119 - t161 * t123 - MDP(16)) * (g(3) * t115 - t159);
t83 = g(1) * t103 + g(2) * t101 + t119 * t148;
t127 = (pkin(3) + t169) * t159;
t126 = -pkin(8) - pkin(7);
t108 = t125 * t151;
t106 = t121 * t151;
t1 = [t137 * MDP(3) + (-g(1) * (-t102 * pkin(4) - t101 * qJ(5)) - g(2) * (t104 * pkin(4) + t103 * qJ(5)) + (g(1) * t126 - g(2) * t131) * t125 + (g(1) * t131 + g(2) * t126) * t121) * MDP(28) + (g(1) * t135 - g(2) * t85) * MDP(34) + (-g(1) * t140 - g(2) * t84) * MDP(35) + t161 * (g(1) * t102 - g(2) * t104) - t160 * (g(1) * t101 - g(2) * t103) + (-t120 * MDP(10) + t115 * MDP(16) + t124 * MDP(9) - t165 * t114 + MDP(2)) * (g(1) * t121 - g(2) * t125); (-g(3) * t124 + t137 * t120) * MDP(9) + (g(3) * t120 + t137 * t124) * MDP(10) + (-g(1) * (-t125 * t153 + t108) - g(2) * (-t121 * t153 + t106) - g(3) * (t116 + t139) + t127) * MDP(28) + t130; (-g(1) * t108 - g(2) * t106 - g(3) * t139 + t127) * MDP(28) + t130; (-g(1) * (-pkin(4) * t103 + qJ(5) * t104) - g(2) * (-pkin(4) * t101 + qJ(5) * t102) - (-pkin(4) * t119 + qJ(5) * t123) * t148) * MDP(28) + t161 * t83 + t160 * (g(1) * t104 + g(2) * t102 + t123 * t148) - t170; -t83 * MDP(28); t170;];
taug  = t1;
