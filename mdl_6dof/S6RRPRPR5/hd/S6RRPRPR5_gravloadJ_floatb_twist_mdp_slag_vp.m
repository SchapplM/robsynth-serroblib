% Calculate Gravitation load on the joints for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:13
% EndTime: 2019-03-09 10:35:17
% DurationCPUTime: 1.07s
% Computational Cost: add. (556->134), mult. (1309->223), div. (0->0), fcn. (1639->14), ass. (0->62)
t140 = sin(pkin(11));
t145 = sin(qJ(2));
t148 = cos(qJ(2));
t184 = cos(pkin(11));
t125 = -t148 * t140 - t145 * t184;
t190 = MDP(19) - MDP(22);
t146 = sin(qJ(1));
t149 = cos(qJ(1));
t154 = -t145 * t140 + t148 * t184;
t143 = cos(pkin(6));
t166 = t125 * t143;
t106 = -t146 * t154 + t149 * t166;
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t141 = sin(pkin(6));
t176 = t141 * t149;
t101 = -t106 * t147 - t144 * t176;
t150 = t143 * t154;
t107 = t146 * t125 + t149 * t150;
t138 = pkin(12) + qJ(6);
t136 = sin(t138);
t137 = cos(t138);
t193 = t101 * t136 + t107 * t137;
t192 = t101 * t137 - t107 * t136;
t170 = t146 * t148;
t173 = t145 * t149;
t122 = -t143 * t170 - t173;
t177 = t141 * t148;
t191 = -g(1) * t122 - g(3) * t177;
t110 = t125 * t149 - t146 * t150;
t117 = t154 * t141;
t189 = -g(1) * t110 - g(2) * t107 - g(3) * t117;
t111 = t146 * t166 + t149 * t154;
t181 = t136 * t147;
t180 = t137 * t147;
t139 = sin(pkin(12));
t179 = t139 * t147;
t178 = t141 * t146;
t142 = cos(pkin(12));
t175 = t142 * t147;
t171 = t146 * t145;
t168 = t148 * t149;
t164 = t143 * t168;
t119 = pkin(2) * t143 * t145 + (-pkin(8) - qJ(3)) * t141;
t135 = pkin(2) * t148 + pkin(1);
t160 = -t119 * t146 + t149 * t135;
t157 = g(1) * t146 - g(2) * t149;
t156 = -t119 * t149 - t146 * t135;
t100 = -t106 * t144 + t147 * t176;
t104 = t111 * t144 - t147 * t178;
t118 = t125 * t141;
t112 = -t118 * t144 - t143 * t147;
t153 = g(1) * t104 + g(2) * t100 + g(3) * t112;
t126 = pkin(2) * t164;
t123 = -t143 * t171 + t168;
t121 = -t143 * t173 - t170;
t120 = -t164 + t171;
t113 = -t118 * t147 + t143 * t144;
t105 = t111 * t147 + t144 * t178;
t99 = t105 * t137 - t110 * t136;
t98 = -t105 * t136 - t110 * t137;
t1 = [t157 * MDP(2) + (-g(1) * t121 - g(2) * t123) * MDP(9) + (-g(1) * t120 - g(2) * t122) * MDP(10) + (-g(1) * t156 - g(2) * t160) * MDP(12) + (g(1) * t101 - g(2) * t105) * MDP(18) + (-g(1) * (-t101 * t142 + t107 * t139) - g(2) * (t105 * t142 - t110 * t139)) * MDP(20) + (-g(1) * (t101 * t139 + t107 * t142) - g(2) * (-t105 * t139 - t110 * t142)) * MDP(21) + (-g(1) * (pkin(3) * t106 - pkin(4) * t101 + t107 * pkin(9) - qJ(5) * t100 + t156) - g(2) * (pkin(3) * t111 + pkin(4) * t105 - pkin(9) * t110 + qJ(5) * t104 + t160)) * MDP(23) + (g(1) * t192 - g(2) * t99) * MDP(29) + (-g(1) * t193 - g(2) * t98) * MDP(30) + t190 * (-g(1) * t100 + g(2) * t104) + (-MDP(11) * t141 + MDP(3)) * (g(1) * t149 + g(2) * t146); (g(2) * t120 + t191) * MDP(9) + (g(3) * t141 * t145 + g(1) * t123 - g(2) * t121) * MDP(10) + (-g(2) * t126 + (g(2) * t171 + t191) * pkin(2)) * MDP(12) + (-g(1) * (t110 * t175 + t111 * t139) - g(2) * (-t106 * t139 + t107 * t175) - g(3) * (t117 * t175 - t118 * t139)) * MDP(20) + (-g(1) * (-t110 * t179 + t111 * t142) - g(2) * (-t106 * t142 - t107 * t179) - g(3) * (-t117 * t179 - t118 * t142)) * MDP(21) + (-g(1) * (t122 * pkin(2) + pkin(9) * t111) - g(2) * (-pkin(2) * t171 - pkin(9) * t106 + t126) - g(3) * (pkin(2) * t177 - pkin(9) * t118) + t189 * (pkin(4) * t147 + qJ(5) * t144 + pkin(3))) * MDP(23) + (-g(1) * (t110 * t180 + t111 * t136) - g(2) * (-t106 * t136 + t107 * t180) - g(3) * (t117 * t180 - t118 * t136)) * MDP(29) + (-g(1) * (-t110 * t181 + t111 * t137) - g(2) * (-t106 * t137 - t107 * t181) - g(3) * (-t117 * t181 - t118 * t137)) * MDP(30) + (t147 * MDP(18) - t190 * t144) * t189; (MDP(12) + MDP(23)) * (-g(3) * t143 - t157 * t141); (-g(1) * (-pkin(4) * t104 + qJ(5) * t105) - g(2) * (-pkin(4) * t100 + qJ(5) * t101) - g(3) * (-pkin(4) * t112 + qJ(5) * t113)) * MDP(23) + t190 * (g(1) * t105 + g(2) * t101 + g(3) * t113) + (MDP(20) * t142 - MDP(21) * t139 + MDP(29) * t137 - MDP(30) * t136 + MDP(18)) * t153; -t153 * MDP(23); (-g(1) * t98 + g(2) * t193 - g(3) * (-t113 * t136 - t117 * t137)) * MDP(29) + (g(1) * t99 + g(2) * t192 - g(3) * (-t113 * t137 + t117 * t136)) * MDP(30);];
taug  = t1;
