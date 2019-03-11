% Calculate Gravitation load on the joints for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:41
% EndTime: 2019-03-09 17:19:42
% DurationCPUTime: 0.78s
% Computational Cost: add. (275->111), mult. (645->154), div. (0->0), fcn. (680->8), ass. (0->57)
t134 = sin(qJ(2));
t133 = sin(qJ(3));
t136 = cos(qJ(5));
t132 = sin(qJ(5));
t137 = cos(qJ(3));
t165 = t132 * t137;
t146 = -t133 * t136 + t165;
t142 = g(3) * t146;
t139 = cos(qJ(1));
t135 = sin(qJ(1));
t138 = cos(qJ(2));
t160 = t135 * t138;
t110 = t133 * t160 + t137 * t139;
t157 = t139 * t133;
t159 = t137 * t138;
t111 = t135 * t159 - t157;
t184 = t110 * t136 - t111 * t132;
t112 = -t135 * t137 + t138 * t157;
t158 = t138 * t139;
t113 = t135 * t133 + t137 * t158;
t96 = t112 * t136 - t113 * t132;
t188 = g(1) * t96 + g(2) * t184 - t134 * t142;
t151 = g(1) * t139 + g(2) * t135;
t180 = t151 * t134;
t101 = -g(3) * t138 + t180;
t145 = t132 * t133 + t136 * t137;
t167 = t110 * t132;
t147 = t111 * t136 + t167;
t170 = g(3) * t134;
t97 = t112 * t132 + t113 * t136;
t186 = (g(1) * t97 + g(2) * t147 + t145 * t170) * MDP(28);
t182 = MDP(16) + MDP(18);
t181 = MDP(17) - MDP(20);
t185 = MDP(10) - MDP(19) + MDP(29);
t176 = pkin(5) * t132;
t175 = g(1) * t135;
t126 = t134 * pkin(8);
t128 = t138 * pkin(2);
t125 = pkin(5) * t136 + pkin(4);
t168 = -pkin(3) - t125;
t164 = t133 * t134;
t163 = t133 * t138;
t162 = t134 * t137;
t161 = t134 * t139;
t156 = -pkin(1) - t128;
t155 = qJ(4) + t176;
t153 = pkin(3) * t159 + qJ(4) * t163 + t126 + t128;
t149 = -t111 * pkin(3) + t139 * pkin(7) - qJ(4) * t110;
t148 = t139 * pkin(1) + pkin(2) * t158 + t113 * pkin(3) + t135 * pkin(7) + pkin(8) * t161;
t95 = g(1) * t112 + g(2) * t110 + g(3) * t164;
t131 = -qJ(6) - pkin(9);
t120 = pkin(8) * t158;
t117 = pkin(8) * t160;
t115 = qJ(4) * t162;
t108 = t112 * pkin(3);
t106 = t110 * pkin(3);
t1 = [t151 * MDP(3) + (-g(1) * t149 - g(2) * (t112 * qJ(4) + t148) - (t156 - t126) * t175) * MDP(21) + (g(1) * t147 - g(2) * t97) * MDP(27) + (g(1) * t184 - g(2) * t96) * MDP(28) + (-g(1) * (-pkin(5) * t167 - t111 * t125 + t149) - g(2) * (t155 * t112 + t113 * t125 + t131 * t161 + t148) - ((-pkin(8) - t131) * t134 + t156) * t175) * MDP(30) + t182 * (g(1) * t111 - g(2) * t113) - t181 * (g(1) * t110 - g(2) * t112) + (t138 * MDP(9) + MDP(2)) * (-g(2) * t139 + t175) - t185 * (-g(2) * t161 + t134 * t175); (-g(1) * t120 - g(2) * t117 - g(3) * t153 + (pkin(3) * t137 + qJ(4) * t133 + pkin(2)) * t180) * MDP(21) + (t138 * t142 - t146 * t180) * MDP(28) + (-g(1) * (t131 * t158 + t120) - g(2) * (t131 * t160 + t117) - g(3) * (t125 * t159 + t163 * t176 + t153) + (-g(3) * t131 + t151 * (t155 * t133 - t168 * t137 + pkin(2))) * t134) * MDP(30) + t185 * (t151 * t138 + t170) + (t145 * MDP(27) - t181 * t133 + t182 * t137 + MDP(9)) * t101; (-g(1) * (qJ(4) * t113 - t108) - g(2) * (qJ(4) * t111 - t106) - g(3) * (-pkin(3) * t164 + t115)) * MDP(21) + t188 * MDP(27) - t186 + (-g(1) * (-t112 * t125 + t155 * t113 - t108) - g(2) * (-t110 * t125 + t155 * t111 - t106) - g(3) * t115 - (pkin(5) * t165 + t168 * t133) * t170) * MDP(30) + t182 * t95 + t181 * (g(1) * t113 + g(2) * t111 + g(3) * t162); -(MDP(21) + MDP(30)) * t95; t186 - (pkin(5) * MDP(30) + MDP(27)) * t188; t101 * MDP(30);];
taug  = t1;
