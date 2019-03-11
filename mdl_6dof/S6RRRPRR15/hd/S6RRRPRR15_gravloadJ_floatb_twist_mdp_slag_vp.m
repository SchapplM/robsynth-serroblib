% Calculate Gravitation load on the joints for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR15_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR15_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR15_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:42
% EndTime: 2019-03-09 20:38:46
% DurationCPUTime: 1.35s
% Computational Cost: add. (656->157), mult. (1802->270), div. (0->0), fcn. (2278->14), ass. (0->70)
t147 = sin(qJ(2));
t148 = sin(qJ(1));
t151 = cos(qJ(2));
t152 = cos(qJ(1));
t181 = cos(pkin(6));
t164 = t152 * t181;
t131 = t148 * t147 - t151 * t164;
t132 = t147 * t164 + t148 * t151;
t146 = sin(qJ(3));
t141 = sin(pkin(7));
t142 = sin(pkin(6));
t183 = cos(qJ(3));
t162 = t141 * t142 * t183;
t143 = cos(pkin(7));
t168 = t143 * t183;
t107 = t131 * t168 + t132 * t146 + t152 * t162;
t175 = t142 * t152;
t120 = -t131 * t141 + t143 * t175;
t145 = sin(qJ(5));
t150 = cos(qJ(5));
t100 = t107 * t145 - t120 * t150;
t174 = t143 * t146;
t108 = -t141 * t146 * t175 - t131 * t174 + t132 * t183;
t144 = sin(qJ(6));
t149 = cos(qJ(6));
t191 = t100 * t144 - t108 * t149;
t190 = t100 * t149 + t108 * t144;
t189 = MDP(16) - MDP(19);
t186 = t107 * t150 + t120 * t145;
t184 = MDP(17) - MDP(20);
t182 = pkin(10) * t141;
t180 = t141 * t145;
t179 = t141 * t150;
t178 = t142 * t147;
t177 = t142 * t148;
t176 = t142 * t151;
t173 = t144 * t145;
t172 = t145 * t149;
t171 = t146 * t147;
t170 = t146 * t151;
t169 = t141 * t178;
t167 = t183 * t147;
t166 = t183 * t151;
t165 = t148 * t181;
t163 = t181 * t141;
t133 = -t152 * t147 - t151 * t165;
t122 = -t133 * t141 + t143 * t177;
t159 = -g(1) * t120 - g(2) * t122;
t134 = -t147 * t165 + t152 * t151;
t111 = -t133 * t168 + t134 * t146 - t148 * t162;
t118 = -t163 * t183 + (-t143 * t166 + t171) * t142;
t157 = g(1) * t111 + g(2) * t107 + g(3) * t118;
t130 = -t141 * t176 + t143 * t181;
t129 = (-t143 * t171 + t166) * t142;
t128 = (t143 * t167 + t170) * t142;
t119 = t146 * t163 + (t143 * t170 + t167) * t142;
t117 = t128 * t145 + t150 * t169;
t116 = t133 * t183 - t134 * t174;
t115 = t133 * t146 + t134 * t168;
t114 = -t131 * t183 - t132 * t174;
t113 = -t131 * t146 + t132 * t168;
t112 = t134 * t183 + (t133 * t143 + t141 * t177) * t146;
t106 = t118 * t145 + t130 * t150;
t104 = t115 * t145 + t134 * t179;
t103 = t113 * t145 + t132 * t179;
t102 = t111 * t145 + t122 * t150;
t101 = t111 * t150 - t122 * t145;
t96 = t102 * t149 + t112 * t144;
t95 = -t102 * t144 + t112 * t149;
t1 = [(g(1) * t148 - g(2) * t152) * MDP(2) + (g(1) * t152 + g(2) * t148) * MDP(3) + (g(1) * t132 - g(2) * t134) * MDP(9) + (-g(1) * t131 - g(2) * t133) * MDP(10) + t159 * MDP(18) + (-g(1) * (-t148 * pkin(1) - t132 * pkin(2) - pkin(3) * t108 + pkin(9) * t175 - qJ(4) * t107) - g(2) * (t152 * pkin(1) + t134 * pkin(2) + t112 * pkin(3) + pkin(9) * t177 + t111 * qJ(4)) + t159 * pkin(10)) * MDP(21) + (g(1) * t100 - g(2) * t102) * MDP(27) + (g(1) * t186 - g(2) * t101) * MDP(28) + (g(1) * t190 - g(2) * t96) * MDP(34) + (-g(1) * t191 - g(2) * t95) * MDP(35) + t184 * (-g(1) * t107 + g(2) * t111) - t189 * (-g(1) * t108 + g(2) * t112); (-g(1) * t133 + g(2) * t131 - g(3) * t176) * MDP(9) + (-g(1) * (t133 * pkin(2) + t116 * pkin(3) + t115 * qJ(4) + t134 * t182) - g(2) * (-t131 * pkin(2) + t114 * pkin(3) + t113 * qJ(4) + t132 * t182) - g(3) * (t129 * pkin(3) + t128 * qJ(4) + (pkin(2) * t151 + t147 * t182) * t142)) * MDP(21) + (-g(1) * t104 - g(2) * t103 - g(3) * t117) * MDP(27) + (-g(1) * (t115 * t150 - t134 * t180) - g(2) * (t113 * t150 - t132 * t180) - g(3) * (t128 * t150 - t145 * t169)) * MDP(28) + (-g(1) * (t104 * t149 + t116 * t144) - g(2) * (t103 * t149 + t114 * t144) - g(3) * (t117 * t149 + t129 * t144)) * MDP(34) + (-g(1) * (-t104 * t144 + t116 * t149) - g(2) * (-t103 * t144 + t114 * t149) - g(3) * (-t117 * t144 + t129 * t149)) * MDP(35) + t184 * (g(1) * t115 + g(2) * t113 + g(3) * t128) - t189 * (g(1) * t116 + g(2) * t114 + g(3) * t129) + (-MDP(18) * t141 + MDP(10)) * (g(1) * t134 + g(2) * t132 + g(3) * t178); (-g(1) * (-t111 * pkin(3) + t112 * qJ(4)) - g(2) * (-t107 * pkin(3) + t108 * qJ(4)) - g(3) * (-t118 * pkin(3) + t119 * qJ(4))) * MDP(21) + (-g(1) * (-t111 * t144 + t112 * t172) - g(2) * (-t107 * t144 + t108 * t172) - g(3) * (-t118 * t144 + t119 * t172)) * MDP(34) + (-g(1) * (-t111 * t149 - t112 * t173) - g(2) * (-t107 * t149 - t108 * t173) - g(3) * (-t118 * t149 - t119 * t173)) * MDP(35) + (-MDP(27) * t145 - MDP(28) * t150 + t184) * (g(1) * t112 + g(2) * t108 + g(3) * t119) + t189 * t157; -t157 * MDP(21); (g(1) * t102 + g(2) * t100 + g(3) * t106) * MDP(28) + (-MDP(34) * t149 + MDP(35) * t144 - MDP(27)) * (g(1) * t101 + g(2) * t186 + g(3) * (t118 * t150 - t130 * t145)); (-g(1) * t95 + g(2) * t191 - g(3) * (-t106 * t144 + t119 * t149)) * MDP(34) + (g(1) * t96 + g(2) * t190 - g(3) * (-t106 * t149 - t119 * t144)) * MDP(35);];
taug  = t1;
