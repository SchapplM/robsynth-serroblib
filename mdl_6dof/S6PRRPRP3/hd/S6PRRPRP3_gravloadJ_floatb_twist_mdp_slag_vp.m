% Calculate Gravitation load on the joints for
% S6PRRPRP3
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
%   see S6PRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:21
% EndTime: 2019-03-08 21:38:23
% DurationCPUTime: 0.77s
% Computational Cost: add. (495->114), mult. (1008->178), div. (0->0), fcn. (1202->12), ass. (0->64)
t178 = MDP(11) - MDP(14) - MDP(24);
t176 = MDP(21) + MDP(23);
t175 = MDP(22) - MDP(25);
t140 = sin(qJ(2));
t142 = cos(qJ(2));
t169 = cos(pkin(10));
t170 = cos(pkin(6));
t153 = t170 * t169;
t168 = sin(pkin(10));
t117 = t140 * t153 + t142 * t168;
t139 = sin(qJ(3));
t141 = cos(qJ(3));
t136 = sin(pkin(6));
t155 = t136 * t169;
t104 = t117 * t139 + t141 * t155;
t152 = t170 * t168;
t119 = -t140 * t152 + t142 * t169;
t154 = t136 * t168;
t106 = t119 * t139 - t141 * t154;
t163 = t136 * t140;
t120 = t139 * t163 - t141 * t170;
t145 = g(1) * t106 + g(2) * t104 + g(3) * t120;
t171 = g(3) * t136;
t134 = pkin(11) + qJ(5);
t132 = sin(t134);
t167 = t132 * t141;
t133 = cos(t134);
t166 = t133 * t141;
t135 = sin(pkin(11));
t165 = t135 * t140;
t164 = t135 * t141;
t162 = t136 * t142;
t137 = cos(pkin(11));
t161 = t137 * t141;
t160 = t141 * t142;
t159 = pkin(2) * t162 + pkin(8) * t163;
t158 = MDP(15) + MDP(26);
t157 = t132 * t162;
t156 = pkin(4) * t135 + pkin(8);
t151 = pkin(3) * t141 + qJ(4) * t139;
t131 = pkin(4) * t137 + pkin(3);
t138 = -pkin(9) - qJ(4);
t150 = t131 * t141 - t138 * t139;
t121 = t139 * t170 + t141 * t163;
t102 = t121 * t132 + t133 * t162;
t105 = t117 * t141 - t139 * t155;
t116 = t140 * t168 - t142 * t153;
t93 = t105 * t132 - t116 * t133;
t107 = t119 * t141 + t139 * t154;
t118 = t140 * t169 + t142 * t152;
t95 = t107 * t132 - t118 * t133;
t148 = g(1) * t95 + g(2) * t93 + g(3) * t102;
t115 = t118 * pkin(2);
t114 = t116 * pkin(2);
t109 = (t132 * t140 + t133 * t160) * t136;
t108 = -t133 * t163 + t141 * t157;
t103 = t121 * t133 - t157;
t101 = -t118 * t166 + t119 * t132;
t100 = -t118 * t167 - t119 * t133;
t99 = -t116 * t166 + t117 * t132;
t98 = -t116 * t167 - t117 * t133;
t96 = t107 * t133 + t118 * t132;
t94 = t105 * t133 + t116 * t132;
t1 = [(-MDP(1) - t158) * g(3); (g(1) * t119 + g(2) * t117 + g(3) * t163) * MDP(4) + (-g(1) * (-t118 * t161 + t119 * t135) - g(2) * (-t116 * t161 + t117 * t135) - (t137 * t160 + t165) * t171) * MDP(12) + (-g(1) * (t118 * t164 + t119 * t137) - g(2) * (t116 * t164 + t117 * t137) - (-t135 * t160 + t137 * t140) * t171) * MDP(13) + (-g(1) * (pkin(8) * t119 - t118 * t151 - t115) - g(2) * (pkin(8) * t117 - t116 * t151 - t114) - g(3) * (t151 * t162 + t159)) * MDP(15) + (-g(1) * (pkin(5) * t101 + qJ(6) * t100 - t118 * t150 + t119 * t156 - t115) - g(2) * (pkin(5) * t99 + qJ(6) * t98 - t116 * t150 + t117 * t156 - t114) - g(3) * (t109 * pkin(5) + t108 * qJ(6) + t159) - (pkin(4) * t165 + t150 * t142) * t171) * MDP(26) + t176 * (-g(1) * t101 - g(2) * t99 - g(3) * t109) + t175 * (g(1) * t100 + g(2) * t98 + g(3) * t108) + (-t141 * MDP(10) + t178 * t139 - MDP(3)) * (-g(1) * t118 - g(2) * t116 + g(3) * t162); (-g(1) * (-pkin(3) * t106 + qJ(4) * t107) - g(2) * (-pkin(3) * t104 + qJ(4) * t105) - g(3) * (-pkin(3) * t120 + qJ(4) * t121)) * MDP(15) + (MDP(26) * t138 + t178) * (g(1) * t107 + g(2) * t105 + g(3) * t121) + (t176 * t133 - t175 * t132 + MDP(26) * (pkin(5) * t133 + qJ(6) * t132 + t131) + MDP(12) * t137 - MDP(13) * t135 + MDP(10)) * t145; -t158 * t145; (-g(1) * (-pkin(5) * t95 + qJ(6) * t96) - g(2) * (-pkin(5) * t93 + qJ(6) * t94) - g(3) * (-pkin(5) * t102 + qJ(6) * t103)) * MDP(26) + t176 * t148 + t175 * (g(1) * t96 + g(2) * t94 + g(3) * t103); -t148 * MDP(26);];
taug  = t1;
