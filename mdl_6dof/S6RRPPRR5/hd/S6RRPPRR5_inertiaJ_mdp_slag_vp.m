% Calculate joint inertia matrix for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR5_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:11:01
% EndTime: 2019-03-09 09:11:04
% DurationCPUTime: 0.72s
% Computational Cost: add. (660->194), mult. (1395->261), div. (0->0), fcn. (1348->8), ass. (0->82)
t138 = sin(pkin(6));
t146 = cos(qJ(2));
t182 = t138 * t146;
t145 = cos(qJ(5));
t141 = sin(qJ(6));
t144 = cos(qJ(6));
t156 = MDP(31) * t144 - MDP(32) * t141;
t153 = MDP(24) + t156;
t191 = t153 * t145;
t142 = sin(qJ(5));
t166 = t142 * MDP(25);
t190 = MDP(16) + MDP(13);
t189 = 0.2e1 * MDP(31);
t188 = 0.2e1 * MDP(32);
t147 = pkin(2) + pkin(3);
t139 = cos(pkin(6));
t187 = pkin(1) * t139;
t143 = sin(qJ(2));
t183 = t138 * t143;
t111 = -t138 * pkin(1) - pkin(2) * t182 - qJ(3) * t183;
t107 = pkin(3) * t182 - t111;
t100 = (pkin(4) * t146 - pkin(9) * t143) * t138 + t107;
t116 = pkin(8) * t182 + t143 * t187;
t128 = t139 * qJ(3);
t110 = t128 + t116;
t106 = -qJ(4) * t182 + t110;
t102 = -pkin(9) * t139 + t106;
t97 = t100 * t145 - t102 * t142;
t95 = -pkin(5) * t182 - t97;
t186 = t141 * t95;
t185 = t144 * t95;
t113 = t139 * t145 + t142 * t183;
t184 = t113 * t142;
t181 = t139 * MDP(8);
t140 = qJ(3) - pkin(9);
t180 = t140 * t141;
t179 = t140 * t144;
t178 = t140 * t145;
t177 = pkin(8) * t183 - t146 * t187;
t176 = MDP(23) * t146;
t175 = MDP(30) * t145;
t114 = -t139 * t142 + t145 * t183;
t104 = t114 * t141 - t144 * t182;
t174 = t104 * MDP(29);
t105 = t114 * t144 + t141 * t182;
t173 = t105 * MDP(26);
t172 = t105 * MDP(28);
t171 = t113 * MDP(30);
t170 = t114 * MDP(20);
t169 = t114 * MDP(21);
t168 = t114 * MDP(25);
t134 = pkin(4) + t147;
t167 = t134 * MDP(24);
t165 = t144 * MDP(26);
t164 = MDP(12) - MDP(17);
t132 = t139 * pkin(2);
t112 = -t132 + t177;
t162 = t144 * t141 * MDP(27);
t161 = -MDP(15) + t166;
t103 = -t139 * pkin(3) - qJ(4) * t183 + t112;
t98 = t100 * t142 + t102 * t145;
t160 = t97 * MDP(24) - t98 * MDP(25);
t159 = -t116 * MDP(10) - t177 * MDP(9);
t158 = -MDP(28) * t144 + MDP(29) * t141;
t117 = pkin(5) * t145 + pkin(10) * t142 + t134;
t108 = t117 * t144 - t141 * t178;
t109 = t117 * t141 + t144 * t178;
t157 = t108 * MDP(31) - t109 * MDP(32);
t155 = -t141 * MDP(31) - t144 * MDP(32);
t101 = t139 * pkin(4) - t103;
t154 = MDP(20) + t158;
t152 = t141 * MDP(28) + t144 * MDP(29) + t155 * pkin(10);
t96 = pkin(10) * t182 + t98;
t99 = pkin(5) * t113 - pkin(10) * t114 + t101;
t93 = -t141 * t96 + t144 * t99;
t94 = t141 * t99 + t144 * t96;
t151 = t93 * MDP(31) - t94 * MDP(32) + t171 + t172 - t174;
t150 = -MDP(22) + t152;
t149 = qJ(3) ^ 2;
t137 = t144 ^ 2;
t135 = t141 ^ 2;
t1 = [t114 ^ 2 * MDP(19) + (t103 ^ 2 + t106 ^ 2 + t107 ^ 2) * MDP(18) + (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) * MDP(14) + MDP(1) + (0.2e1 * MDP(6) * t183 + t181) * t139 + (-0.2e1 * t104 * MDP(27) + t173) * t105 + 0.2e1 * (MDP(7) * t139 + t169) * t182 + (-0.2e1 * MDP(22) * t182 - 0.2e1 * t170 + t171 + 0.2e1 * t172 - 0.2e1 * t174) * t113 + (t104 * t95 + t113 * t93) * t189 + (t105 * t95 - t113 * t94) * t188 + 0.2e1 * (t113 * MDP(24) + t168) * t101 + 0.2e1 * (-t112 * MDP(11) + t110 * MDP(13) - t103 * MDP(15) + t106 * MDP(16) + t159) * t139 + 0.2e1 * ((t112 * MDP(12) - t111 * MDP(13) + t107 * MDP(16) - t103 * MDP(17)) * t143 + (-t111 * MDP(11) + t110 * MDP(12) + t107 * MDP(15) - t106 * MDP(17) + t160) * t146) * t138 + (t143 ^ 2 * MDP(4) + (0.2e1 * MDP(5) * t143 + t176) * t146 + 0.2e1 * (-MDP(10) * t143 + MDP(9) * t146) * pkin(1)) * t138 ^ 2; t134 * t168 + (t139 * t147 - t103) * MDP(15) + (qJ(3) * t106 - t103 * t147) * MDP(18) + t181 + (0.2e1 * t132 - t177) * MDP(11) + (-pkin(2) * t112 + qJ(3) * t110) * MDP(14) + (t157 + t167) * t113 + (t101 * MDP(24) + t151 - t170) * t145 + ((t104 * t144 + t105 * t141) * MDP(27) + (t105 * t140 - t185) * MDP(32) + (t104 * t140 - t186) * MDP(31) - t101 * MDP(25) - t105 * t165 - t114 * MDP(19) + t154 * t113) * t142 + ((-pkin(2) * MDP(12) + t147 * MDP(17) + MDP(6)) * t143 + (-qJ(4) * MDP(16) - t142 * MDP(21) - t145 * MDP(22) + MDP(7) + (-MDP(24) * t142 - MDP(25) * t145) * t140 + t164 * qJ(3)) * t146) * t138 + t159 + t190 * (0.2e1 * t128 + t116); MDP(8) + 0.2e1 * pkin(2) * MDP(11) + (pkin(2) ^ 2 + t149) * MDP(14) + (t147 ^ 2 + t149) * MDP(18) + 0.2e1 * t147 * MDP(15) - 0.2e1 * t134 * t166 + 0.2e1 * t190 * qJ(3) + (t137 * MDP(26) - t179 * t188 - t180 * t189 + MDP(19) - 0.2e1 * t162) * t142 ^ 2 + (t108 * t189 - t109 * t188 + 0.2e1 * t154 * t142 + 0.2e1 * t167 + t175) * t145; t112 * MDP(14) + t103 * MDP(18) - t168 + (-MDP(15) - MDP(11)) * t139 + t164 * t183 - t153 * t113; -pkin(2) * MDP(14) - MDP(18) * t147 - MDP(11) + t161 - t191; MDP(14) + MDP(18); t107 * MDP(18) + (-t104 * t145 - t141 * t184) * MDP(31) + (-t105 * t145 - t144 * t184) * MDP(32) + (t143 * MDP(16) + (t145 * MDP(24) - t161) * t146) * t138; 0; 0; MDP(18); t169 + t138 * t176 + t141 * t173 + (-t104 * t141 + t105 * t144) * MDP(27) + (-pkin(5) * t104 - t185) * MDP(31) + (-pkin(5) * t105 + t186) * MDP(32) + t150 * t113 + t160; (-t140 * MDP(25) + t150) * t145 + (-MDP(21) - t140 * MDP(24) - t141 * t165 + (t135 - t137) * MDP(27) + (pkin(5) * t141 - t179) * MDP(31) + (pkin(5) * t144 + t180) * MDP(32)) * t142; 0; -t166 + t191; MDP(26) * t135 + 0.2e1 * pkin(5) * t156 + MDP(23) + 0.2e1 * t162; t151; t158 * t142 + t157 + t175; -t156; t155 * t142; t152; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
