% Calculate joint inertia matrix for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR8_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:52
% EndTime: 2019-03-09 05:24:54
% DurationCPUTime: 0.63s
% Computational Cost: add. (820->172), mult. (1567->253), div. (0->0), fcn. (1651->8), ass. (0->78)
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t182 = -qJ(5) - pkin(8);
t138 = t182 * t152;
t139 = t182 * t155;
t149 = sin(pkin(10));
t150 = cos(pkin(10));
t117 = t150 * t138 + t139 * t149;
t132 = t149 * t155 + t150 * t152;
t105 = -pkin(9) * t132 + t117;
t118 = t149 * t138 - t150 * t139;
t131 = -t149 * t152 + t150 * t155;
t106 = pkin(9) * t131 + t118;
t151 = sin(qJ(6));
t154 = cos(qJ(6));
t112 = -t154 * t131 + t132 * t151;
t113 = t131 * t151 + t132 * t154;
t164 = t113 * MDP(25) - t112 * MDP(26) + (t105 * t154 - t106 * t151) * MDP(28) - (t105 * t151 + t106 * t154) * MDP(29);
t188 = MDP(19) * t152 + MDP(20) * t155;
t190 = t152 * MDP(16) + t155 * MDP(17) - pkin(8) * t188 + t164;
t156 = cos(qJ(3));
t124 = t132 * t156;
t177 = t155 * t156;
t179 = t152 * t156;
t126 = -t149 * t179 + t150 * t177;
t102 = t154 * t124 + t126 * t151;
t104 = -t124 * t151 + t126 * t154;
t189 = t104 * MDP(25) - t102 * MDP(26);
t187 = 2 * MDP(21);
t186 = -2 * MDP(24);
t185 = 0.2e1 * MDP(29);
t184 = (pkin(1) * MDP(6));
t183 = pkin(4) * t149;
t153 = sin(qJ(3));
t123 = t132 * t153;
t125 = t131 * t153;
t101 = -t123 * t154 - t125 * t151;
t103 = -t123 * t151 + t125 * t154;
t181 = t101 * MDP(28) - t103 * MDP(29);
t180 = t152 * t155;
t157 = -pkin(1) - pkin(7);
t178 = t152 * t157;
t176 = t155 * t157;
t137 = pkin(3) * t153 - pkin(8) * t156 + qJ(2);
t130 = t155 * t137;
t111 = -qJ(5) * t177 + t130 + (pkin(4) - t178) * t153;
t166 = t153 * t176;
t116 = t166 + (-qJ(5) * t156 + t137) * t152;
t96 = t149 * t111 + t150 * t116;
t172 = MDP(23) * t113;
t171 = MDP(28) * t112;
t142 = pkin(4) * t150 + pkin(5);
t127 = t142 * t154 - t151 * t183;
t170 = MDP(28) * t127;
t128 = t142 * t151 + t154 * t183;
t169 = t128 * MDP(29);
t168 = MDP(18) + MDP(27);
t167 = t153 * MDP(27) + t189;
t143 = -pkin(4) * t155 - pkin(3);
t165 = MDP(15) * t180;
t95 = t150 * t111 - t116 * t149;
t89 = pkin(5) * t153 - pkin(9) * t126 + t95;
t92 = -pkin(9) * t124 + t96;
t86 = -t151 * t92 + t154 * t89;
t133 = pkin(4) * t179 - t156 * t157;
t87 = t151 * t89 + t154 * t92;
t163 = MDP(16) * t155 - MDP(17) * t152;
t162 = MDP(19) * t155 - MDP(20) * t152;
t160 = MDP(22) * t143 + MDP(29) * t113 + t171;
t148 = t156 ^ 2;
t147 = t155 ^ 2;
t146 = t153 ^ 2;
t145 = t152 ^ 2;
t122 = -pkin(5) * t131 + t143;
t120 = t137 * t152 + t166;
t119 = -t153 * t178 + t130;
t115 = pkin(5) * t124 + t133;
t1 = [(t133 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(22) + MDP(1) + t168 * t146 + (MDP(23) * t104 + t102 * t186) * t104 + ((-2 * MDP(4) + t184) * pkin(1)) + (0.2e1 * t156 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t147 * MDP(14) + MDP(7) - 0.2e1 * t165) * t148 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t163) * t156 + t189) * t153 + 0.2e1 * (t119 * t153 - t148 * t178) * MDP(19) + 0.2e1 * (-t120 * t153 - t148 * t176) * MDP(20) + (-t124 * t96 - t95 * t126) * t187 + 0.2e1 * (t102 * t115 + t153 * t86) * MDP(28) + (t104 * t115 - t153 * t87) * t185; MDP(4) - t184 + (t123 * t126 - t124 * t125) * MDP(21) + (-t123 * t95 + t125 * t96 - t133 * t156) * MDP(22) + (t101 * t153 - t102 * t156) * MDP(28) + (-t103 * t153 - t104 * t156) * MDP(29) + t188 * (-t146 - t148); MDP(6) + (t123 ^ 2 + t125 ^ 2 + t148) * MDP(22); (-t117 * t126 - t118 * t124 + t131 * t96 - t132 * t95) * MDP(21) + (t117 * t95 + t118 * t96 + t133 * t143) * MDP(22) + t104 * t172 + (-t102 * t113 - t104 * t112) * MDP(24) + (t102 * t122 + t112 * t115) * MDP(28) + (t104 * t122 + t113 * t115) * MDP(29) + (-t157 * MDP(13) - MDP(10) + t190) * t153 + (MDP(9) + t157 * MDP(12) + MDP(14) * t180 + (-t145 + t147) * MDP(15) + (-pkin(3) * t152 + t176) * MDP(19) + (-pkin(3) * t155 - t178) * MDP(20)) * t156; -t153 * MDP(13) + (t123 * t132 + t125 * t131) * MDP(21) + (-t117 * t123 + t118 * t125) * MDP(22) + (MDP(12) - t160 + t162) * t156; MDP(11) + t145 * MDP(14) + 0.2e1 * t165 + (-t117 * t132 + t118 * t131) * t187 + (t117 ^ 2 + t118 ^ 2 + t143 ^ 2) * MDP(22) + 0.2e1 * t122 * t171 + 0.2e1 * t162 * pkin(3) + (t112 * t186 + t122 * t185 + t172) * t113; t153 * MDP(18) + t119 * MDP(19) - t120 * MDP(20) + (t127 * t153 + t86) * MDP(28) + (-t128 * t153 - t87) * MDP(29) + t163 * t156 + ((-t124 * t149 - t126 * t150) * MDP(21) + (t149 * t96 + t150 * t95) * MDP(22)) * pkin(4) + t167; -t188 * t153 + (-t123 * t150 + t125 * t149) * MDP(22) * pkin(4) + t181; ((t131 * t149 - t132 * t150) * MDP(21) + (t117 * t150 + t118 * t149) * MDP(22)) * pkin(4) + t190; (t149 ^ 2 + t150 ^ 2) * MDP(22) * pkin(4) ^ 2 + 0.2e1 * t170 - 0.2e1 * t169 + t168; MDP(22) * t133 + t102 * MDP(28) + t104 * MDP(29); -t156 * MDP(22); t160; 0; MDP(22); t86 * MDP(28) - t87 * MDP(29) + t167; t181; t164; MDP(27) - t169 + t170; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
