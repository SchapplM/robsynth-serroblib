% Calculate Gravitation load on the joints for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:50
% EndTime: 2021-01-16 02:38:54
% DurationCPUTime: 0.78s
% Computational Cost: add. (440->103), mult. (789->160), div. (0->0), fcn. (911->12), ass. (0->58)
t128 = sin(qJ(3));
t131 = cos(qJ(3));
t124 = sin(pkin(6));
t129 = sin(qJ(2));
t159 = t124 * t129;
t168 = cos(pkin(6));
t177 = -t128 * t159 + t168 * t131;
t123 = sin(pkin(10));
t148 = t123 * t168;
t167 = cos(pkin(10));
t170 = cos(qJ(2));
t105 = t129 * t148 - t167 * t170;
t158 = t124 * t131;
t176 = t105 * t128 + t123 * t158;
t122 = qJ(3) + pkin(11);
t120 = sin(t122);
t121 = cos(t122);
t102 = t120 * t159 - t168 * t121;
t143 = t168 * t167;
t107 = t123 * t170 + t129 * t143;
t147 = t124 * t167;
t96 = t107 * t120 + t121 * t147;
t160 = t123 * t124;
t98 = t105 * t120 + t121 * t160;
t138 = g(1) * t98 - g(2) * t96 - g(3) * t102;
t130 = cos(qJ(5));
t118 = pkin(5) * t130 + pkin(4);
t125 = -qJ(6) - pkin(9);
t142 = -t118 * t121 + t120 * t125;
t175 = -MDP(13) + MDP(25);
t174 = MDP(21) + MDP(23);
t173 = MDP(22) + MDP(24);
t169 = g(3) * t124;
t127 = sin(qJ(5));
t166 = t105 * t127;
t162 = t121 * t127;
t161 = t121 * t130;
t157 = t127 * t129;
t108 = t167 * t129 + t170 * t148;
t119 = pkin(3) * t131 + pkin(2);
t126 = qJ(4) + pkin(8);
t156 = -t105 * t126 - t108 * t119;
t155 = MDP(15) + MDP(26);
t152 = t124 * t170;
t151 = t127 * t170;
t150 = t130 * t170;
t149 = g(3) * (t119 * t152 + t126 * t159);
t95 = t105 * t121 - t120 * t160;
t145 = t176 * pkin(3);
t141 = t177 * pkin(3);
t136 = -t107 * t128 - t131 * t147;
t133 = t136 * pkin(3);
t106 = t123 * t129 - t170 * t143;
t132 = g(1) * t108 + g(2) * t106 - g(3) * t152;
t103 = t168 * t120 + t121 * t159;
t100 = t106 * t119;
t97 = t107 * t121 - t120 * t147;
t1 = [(-MDP(1) - t155) * g(3); (-g(1) * t156 - g(2) * (t107 * t126 - t100) - t149) * MDP(15) + (-g(1) * (-pkin(5) * t166 + t142 * t108 + t156) - g(2) * (-t100 + (pkin(5) * t127 + t126) * t107 + t142 * t106) - t149 - (pkin(5) * t157 - t142 * t170) * t169) * MDP(26) + t174 * (-g(1) * (-t108 * t161 - t166) - g(2) * (-t106 * t161 + t107 * t127) - (t121 * t150 + t157) * t169) + t173 * (-g(1) * (-t105 * t130 + t108 * t162) - g(2) * (t106 * t162 + t107 * t130) - (-t121 * t151 + t129 * t130) * t169) + (MDP(4) - MDP(14)) * (-g(1) * t105 + g(2) * t107 + g(3) * t159) + (t131 * MDP(10) - t128 * MDP(11) + t121 * MDP(12) + t175 * t120 + MDP(3)) * t132; (-g(1) * t176 - g(2) * t136 - g(3) * t177) * MDP(10) + (-g(1) * (t105 * t131 - t128 * t160) - g(2) * (-t107 * t131 + t128 * t147) - g(3) * (-t168 * t128 - t129 * t158)) * MDP(11) + (-g(1) * t145 - g(2) * t133 - g(3) * t141) * MDP(15) + (-g(1) * (t118 * t98 + t125 * t95 + t145) - g(2) * (-t96 * t118 - t97 * t125 + t133) - g(3) * (-t102 * t118 - t103 * t125 + t141)) * MDP(26) + t175 * (g(1) * t95 - g(2) * t97 - g(3) * t103) + (t173 * t127 - t174 * t130 - MDP(12)) * t138; -t155 * t132; t173 * (-g(1) * (-t108 * t127 + t130 * t95) - g(2) * (-t106 * t127 - t130 * t97) - g(3) * (-t103 * t130 + t124 * t151)) + (MDP(26) * pkin(5) + t174) * (-g(1) * (t108 * t130 + t127 * t95) - g(2) * (t106 * t130 - t127 * t97) - g(3) * (-t103 * t127 - t124 * t150)); t138 * MDP(26);];
taug = t1;
