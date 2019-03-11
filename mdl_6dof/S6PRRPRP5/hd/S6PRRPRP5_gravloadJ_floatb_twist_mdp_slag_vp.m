% Calculate Gravitation load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:17
% EndTime: 2019-03-08 21:49:19
% DurationCPUTime: 0.66s
% Computational Cost: add. (382->104), mult. (976->160), div. (0->0), fcn. (1163->10), ass. (0->61)
t189 = -MDP(10) + MDP(13) - MDP(24);
t148 = sin(qJ(3));
t188 = -qJ(4) * t148 - pkin(2);
t186 = MDP(11) - MDP(14);
t167 = MDP(21) + MDP(23);
t185 = MDP(22) - MDP(25);
t184 = pkin(4) + pkin(8);
t180 = cos(pkin(6));
t179 = cos(pkin(10));
t145 = sin(pkin(10));
t149 = sin(qJ(2));
t152 = cos(qJ(2));
t160 = t180 * t179;
t129 = t145 * t149 - t152 * t160;
t151 = cos(qJ(3));
t177 = t129 * t151;
t163 = t145 * t180;
t131 = t149 * t179 + t152 * t163;
t176 = t131 * t151;
t146 = sin(pkin(6));
t175 = t146 * t149;
t174 = t146 * t151;
t173 = t146 * t152;
t147 = sin(qJ(5));
t172 = t147 * t148;
t150 = cos(qJ(5));
t171 = t148 * t150;
t170 = t148 * t152;
t169 = t151 * t152;
t168 = MDP(15) + MDP(26);
t166 = t150 * t173;
t165 = -pkin(3) * t177 + t129 * t188;
t164 = -pkin(3) * t176 + t131 * t188;
t162 = t146 * t179;
t161 = pkin(2) * t173 + pkin(8) * t175 + (pkin(3) * t169 + qJ(4) * t170) * t146;
t132 = -t149 * t163 + t152 * t179;
t114 = t132 * t148 - t145 * t174;
t100 = -t114 * t150 + t131 * t147;
t133 = t148 * t175 - t151 * t180;
t116 = t133 * t150 + t147 * t173;
t130 = t145 * t152 + t149 * t160;
t112 = t130 * t148 + t151 * t162;
t98 = -t112 * t150 + t129 * t147;
t158 = g(1) * t100 + g(2) * t98 - g(3) * t116;
t97 = g(1) * t114 + g(2) * t112 + g(3) * t133;
t134 = t148 * t180 + t149 * t174;
t128 = t133 * pkin(3);
t121 = (t147 * t170 + t149 * t150) * t146;
t120 = t147 * t175 - t148 * t166;
t117 = -t133 * t147 + t166;
t115 = t145 * t146 * t148 + t132 * t151;
t113 = t130 * t151 - t148 * t162;
t111 = t114 * pkin(3);
t110 = t112 * pkin(3);
t107 = -t131 * t172 + t132 * t150;
t106 = t131 * t171 + t132 * t147;
t105 = -t129 * t172 + t130 * t150;
t104 = t129 * t171 + t130 * t147;
t101 = t114 * t147 + t131 * t150;
t99 = t112 * t147 + t129 * t150;
t1 = [(-MDP(1) - t168) * g(3); (-g(1) * (pkin(8) * t132 + t164) - g(2) * (pkin(8) * t130 + t165) - g(3) * t161) * MDP(15) + (-g(1) * (pkin(5) * t107 - pkin(9) * t176 + qJ(6) * t106 + t132 * t184 + t164) - g(2) * (pkin(5) * t105 - pkin(9) * t177 + qJ(6) * t104 + t184 * t130 + t165) - g(3) * (t121 * pkin(5) + t120 * qJ(6) + (pkin(4) * t149 + pkin(9) * t169) * t146 + t161)) * MDP(26) + t167 * (-g(1) * t107 - g(2) * t105 - g(3) * t121) + t185 * (g(1) * t106 + g(2) * t104 + g(3) * t120) + (MDP(4) - MDP(12)) * (g(1) * t132 + g(2) * t130 + g(3) * t175) + (t186 * t148 + t189 * t151 - MDP(3)) * (-g(1) * t131 - g(2) * t129 + g(3) * t173); (-g(1) * (qJ(4) * t115 - t111) - g(2) * (qJ(4) * t113 - t110) - g(3) * (qJ(4) * t134 - t128)) * MDP(15) + (-g(1) * (-pkin(9) * t114 - t111) - g(2) * (-pkin(9) * t112 - t110) - g(3) * (-pkin(9) * t133 - t128)) * MDP(26) - t189 * t97 + (-t185 * t150 - t167 * t147 + (-pkin(5) * t147 + qJ(6) * t150 - qJ(4)) * MDP(26) + t186) * (g(1) * t115 + g(2) * t113 + g(3) * t134); -t168 * t97; (-g(1) * (-pkin(5) * t100 + qJ(6) * t101) - g(2) * (-pkin(5) * t98 + qJ(6) * t99) - g(3) * (pkin(5) * t116 - qJ(6) * t117)) * MDP(26) + t167 * t158 + t185 * (g(1) * t101 + g(2) * t99 - g(3) * t117); -t158 * MDP(26);];
taug  = t1;
