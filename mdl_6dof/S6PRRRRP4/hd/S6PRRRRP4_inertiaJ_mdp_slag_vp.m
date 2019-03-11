% Calculate joint inertia matrix for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:49
% EndTime: 2019-03-09 00:17:52
% DurationCPUTime: 1.02s
% Computational Cost: add. (917->207), mult. (1881->296), div. (0->0), fcn. (1990->10), ass. (0->79)
t147 = sin(qJ(4));
t151 = cos(qJ(4));
t186 = -pkin(10) - pkin(9);
t130 = t186 * t151;
t146 = sin(qJ(5));
t150 = cos(qJ(5));
t165 = t186 * t147;
t112 = -t146 * t130 - t150 * t165;
t113 = -t150 * t130 + t146 * t165;
t176 = t150 * t151;
t125 = t146 * t147 - t176;
t126 = t146 * t151 + t150 * t147;
t167 = MDP(25) - MDP(28);
t168 = MDP(24) + MDP(26);
t154 = t126 * MDP(21) - t125 * MDP(22) - t168 * t112 - t167 * t113;
t198 = t154 - (t147 * MDP(17) + t151 * MDP(18)) * pkin(9) + t147 * MDP(14) + t151 * MDP(15);
t170 = t146 * MDP(25);
t195 = (t150 * MDP(24) - t170) * pkin(4);
t148 = sin(qJ(3));
t118 = t126 * t148;
t179 = t147 * t148;
t119 = -t146 * t179 + t148 * t176;
t173 = t119 * MDP(21) - t118 * MDP(22);
t152 = cos(qJ(3));
t129 = -t152 * pkin(3) - t148 * pkin(9) - pkin(2);
t124 = t151 * t129;
t177 = t148 * t151;
t185 = pkin(8) * t147;
t103 = -pkin(10) * t177 + t124 + (-pkin(4) - t185) * t152;
t183 = pkin(8) * t152;
t166 = t151 * t183;
t107 = t166 + (-pkin(10) * t148 + t129) * t147;
t93 = t150 * t103 - t146 * t107;
t193 = t93 * MDP(24) + t173;
t191 = -2 * MDP(20);
t190 = 0.2e1 * MDP(25);
t189 = 2 * MDP(26);
t188 = 2 * MDP(27);
t187 = 2 * MDP(28);
t184 = pkin(8) * t151;
t182 = t152 * pkin(5);
t94 = t146 * t103 + t150 * t107;
t144 = sin(pkin(6));
t149 = sin(qJ(2));
t181 = t144 * t149;
t153 = cos(qJ(2));
t180 = t144 * t153;
t178 = t147 * t151;
t175 = t152 * qJ(6);
t128 = pkin(4) * t179 + t148 * pkin(8);
t172 = MDP(11) * t148;
t171 = t119 * MDP(19);
t169 = MDP(16) + MDP(23);
t137 = -t151 * pkin(4) - pkin(3);
t164 = MDP(13) * t178;
t145 = cos(pkin(6));
t121 = t145 * t148 + t152 * t181;
t105 = -t121 * t147 - t151 * t180;
t106 = t121 * t151 - t147 * t180;
t95 = -t150 * t105 + t146 * t106;
t96 = t146 * t105 + t150 * t106;
t163 = -t167 * t96 - t168 * t95;
t161 = pkin(5) * t189 + MDP(23);
t159 = t151 * MDP(14) - t147 * MDP(15);
t157 = t151 * MDP(17) - t147 * MDP(18);
t142 = t151 ^ 2;
t141 = t148 ^ 2;
t140 = t147 ^ 2;
t138 = t146 * pkin(4);
t135 = t150 * pkin(4) + pkin(5);
t133 = t138 + qJ(6);
t120 = -t145 * t152 + t148 * t181;
t115 = t147 * t129 + t166;
t114 = -t147 * t183 + t124;
t100 = t125 * pkin(5) - t126 * qJ(6) + t137;
t97 = t118 * pkin(5) - t119 * qJ(6) + t128;
t88 = -t93 + t182;
t87 = -t175 + t94;
t1 = [MDP(1) + (t120 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(29); (-t105 * t152 + t120 * t179) * MDP(17) + (t106 * t152 + t120 * t177) * MDP(18) + (-t96 * t118 + t95 * t119) * MDP(27) + (t120 * t97 + t96 * t87 + t95 * t88) * MDP(29) + (-t149 * MDP(4) + (MDP(10) * t152 + MDP(3) - t172) * t153) * t144 + t168 * (t120 * t118 + t95 * t152) + t167 * (t120 * t119 + t96 * t152); MDP(2) - 0.2e1 * pkin(2) * t172 + (t87 ^ 2 + t88 ^ 2 + t97 ^ 2) * MDP(29) + t169 * t152 ^ 2 + (t118 * t191 + t171) * t119 + (t142 * MDP(12) + MDP(5) - 0.2e1 * t164) * t141 + 0.2e1 * (pkin(2) * MDP(10) + (MDP(6) - t159) * t148 - t173) * t152 + 0.2e1 * (-t114 * t152 + t141 * t185) * MDP(17) + 0.2e1 * (t115 * t152 + t141 * t184) * MDP(18) + 0.2e1 * (t128 * t118 - t93 * t152) * MDP(24) + (t128 * t119 + t94 * t152) * t190 + (t97 * t118 + t88 * t152) * t189 + (-t87 * t118 + t88 * t119) * t188 + (-t97 * t119 - t87 * t152) * t187; -t121 * MDP(11) + (-t96 * t125 + t95 * t126) * MDP(27) + (t95 * t112 + t96 * t113) * MDP(29) + (t100 * MDP(29) + t168 * t125 + t167 * t126 - MDP(10) - t157) * t120; t126 * t171 + (-t126 * t118 - t119 * t125) * MDP(20) + (t137 * t118 + t128 * t125) * MDP(24) + (t137 * t119 + t128 * t126) * MDP(25) + (t100 * t118 + t97 * t125) * MDP(26) + (t112 * t119 - t113 * t118 - t87 * t125 + t88 * t126) * MDP(27) + (-t100 * t119 - t97 * t126) * MDP(28) + (t97 * t100 + t88 * t112 + t87 * t113) * MDP(29) + (-pkin(8) * MDP(11) + MDP(8) - t198) * t152 + (MDP(7) - pkin(8) * MDP(10) + MDP(12) * t178 + (-t140 + t142) * MDP(13) + (-pkin(3) * t147 - t184) * MDP(17) + (-pkin(3) * t151 + t185) * MDP(18)) * t148; MDP(9) + t140 * MDP(12) + 0.2e1 * t164 + (t100 ^ 2 + t112 ^ 2 + t113 ^ 2) * MDP(29) + (MDP(19) * t126 - 0.2e1 * t100 * MDP(28) + t112 * t188 + t125 * t191 + t137 * t190) * t126 + 0.2e1 * t157 * pkin(3) + 0.2e1 * (t137 * MDP(24) + t100 * MDP(26) - MDP(27) * t113) * t125; t105 * MDP(17) - t106 * MDP(18) + (t96 * t133 - t95 * t135) * MDP(29) + t163; t114 * MDP(17) - t115 * MDP(18) + t93 * MDP(26) + (-t133 * t118 - t135 * t119) * MDP(27) + (t87 * t133 - t88 * t135) * MDP(29) + t159 * t148 + ((-pkin(5) - t135) * MDP(26) + (-qJ(6) - t133) * MDP(28) - t195 - t169) * t152 - t167 * t94 + t193; (-t133 * t125 - t135 * t126) * MDP(27) + (-t112 * t135 + t113 * t133) * MDP(29) + t198; (t133 ^ 2 + t135 ^ 2) * MDP(29) + 0.2e1 * t195 + t135 * t189 + t133 * t187 + t169; (-t95 * pkin(5) + t96 * qJ(6)) * MDP(29) + t163; -t152 * MDP(23) - t94 * MDP(25) + (t93 - 0.2e1 * t182) * MDP(26) + (-pkin(5) * t119 - t118 * qJ(6)) * MDP(27) + (-0.2e1 * t175 + t94) * MDP(28) + (-t88 * pkin(5) + t87 * qJ(6)) * MDP(29) + t193; (-pkin(5) * t126 - t125 * qJ(6)) * MDP(27) + (-t112 * pkin(5) + t113 * qJ(6)) * MDP(29) + t154; (0.2e1 * qJ(6) + t138) * MDP(28) + (t135 * pkin(5) + t133 * qJ(6)) * MDP(29) + (t168 * t150 - t170) * pkin(4) + t161; qJ(6) * t187 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29) + t161; t95 * MDP(29); t152 * MDP(26) + t119 * MDP(27) + t88 * MDP(29); t126 * MDP(27) + t112 * MDP(29); -t135 * MDP(29) - MDP(26); -MDP(29) * pkin(5) - MDP(26); MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
