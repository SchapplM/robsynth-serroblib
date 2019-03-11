% Calculate joint inertia matrix for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP4_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:41
% EndTime: 2019-03-09 10:01:42
% DurationCPUTime: 0.62s
% Computational Cost: add. (953->181), mult. (1554->247), div. (0->0), fcn. (1498->6), ass. (0->75)
t154 = sin(pkin(9));
t155 = cos(pkin(9));
t156 = sin(qJ(4));
t158 = cos(qJ(4));
t128 = -t154 * t156 + t155 * t158;
t165 = t154 * t158 + t155 * t156;
t209 = t128 ^ 2 + t165 ^ 2;
t196 = pkin(2) + pkin(8);
t181 = -qJ(5) - t196;
t133 = t181 * t156;
t171 = t181 * t158;
t119 = t133 * t154 - t155 * t171;
t121 = t155 * t133 + t154 * t171;
t208 = -t119 * t128 + t121 * t165;
t207 = MDP(14) * pkin(7) + MDP(11);
t206 = MDP(22) + MDP(25);
t205 = MDP(23) + MDP(27);
t204 = -t119 * MDP(24) + t121 * MDP(26);
t159 = cos(qJ(2));
t202 = 0.2e1 * t159;
t201 = 2 * MDP(20);
t200 = 2 * MDP(21);
t199 = 2 * MDP(22);
t198 = 0.2e1 * MDP(24);
t197 = 2 * MDP(25);
t148 = t159 * pkin(7);
t193 = pkin(2) * MDP(14);
t157 = sin(qJ(2));
t135 = (pkin(3) + pkin(7)) * t157;
t190 = t135 * t156;
t189 = t156 * t159;
t188 = t158 * t159;
t187 = t159 * qJ(3);
t145 = pkin(4) * t156 + qJ(3);
t132 = t158 * t135;
t177 = -qJ(3) * t157 - pkin(1);
t126 = -t159 * t196 + t177;
t176 = qJ(5) * t159 - t126;
t111 = pkin(4) * t157 + t156 * t176 + t132;
t113 = -t158 * t176 + t190;
t103 = t111 * t154 + t113 * t155;
t136 = pkin(3) * t159 + t148;
t186 = MDP(14) * pkin(7) ^ 2;
t185 = MDP(24) * t165;
t184 = MDP(26) * t128;
t114 = pkin(5) * t165 - qJ(6) * t128 + t145;
t183 = MDP(27) * t114;
t139 = pkin(4) * t154 + qJ(6);
t182 = t139 * MDP(26);
t180 = t119 ^ 2 + t121 ^ 2;
t100 = qJ(6) * t157 + t103;
t124 = pkin(4) * t188 + t136;
t178 = t158 * t156 * MDP(16);
t122 = t154 * t189 - t155 * t188;
t123 = t165 * t159;
t175 = -t119 * t123 + t121 * t122;
t102 = t111 * t155 - t113 * t154;
t101 = -pkin(5) * t157 - t102;
t169 = t100 * t165 - t101 * t128;
t168 = t102 * t128 + t103 * t165;
t143 = pkin(4) * t155 + pkin(5);
t167 = t128 * t143 + t139 * t165;
t166 = t128 * t155 + t154 * t165;
t164 = -MDP(17) * t156 - MDP(18) * t158;
t163 = -t157 * t196 + t187;
t162 = t158 * MDP(20) - t156 * MDP(21) + t128 * MDP(24) + MDP(26) * t165;
t153 = t159 ^ 2;
t152 = t158 ^ 2;
t151 = t157 ^ 2;
t150 = t156 ^ 2;
t134 = -pkin(2) * t159 + t177;
t117 = t126 * t158 + t190;
t116 = -t126 * t156 + t132;
t104 = -pkin(5) * t122 + qJ(6) * t123 + t124;
t1 = [MDP(1) + pkin(1) * MDP(9) * t202 + (t102 ^ 2 + t103 ^ 2 + t124 ^ 2) * MDP(23) + (t100 ^ 2 + t101 ^ 2 + t104 ^ 2) * MDP(27) + (MDP(12) * t202 + MDP(14) * t134) * t134 + (MDP(15) * t150 + 0.2e1 * t178 + t186) * t153 + (MDP(19) + MDP(4) + t186) * t151 + 0.2e1 * (-pkin(1) * MDP(10) - t134 * MDP(13) + (MDP(5) + t164) * t159) * t157 + (t116 * t157 + t136 * t188) * t201 + (-t117 * t157 - t136 * t189) * t200 + (t102 * t123 + t103 * t122) * t199 + (-t101 * t157 - t104 * t122) * t198 + (t100 * t122 - t101 * t123) * t197 + 0.2e1 * (t100 * t157 + t104 * t123) * MDP(26) + 0.2e1 * (t151 + t153) * MDP(11) * pkin(7); -t156 * MDP(15) * t188 + (t136 * t156 + t158 * t163) * MDP(20) + (t136 * t158 - t156 * t163) * MDP(21) + (-t168 + t175) * MDP(22) + (-t102 * t119 + t103 * t121 + t124 * t145) * MDP(23) + (t104 * t165 - t114 * t122) * MDP(24) + (-t169 + t175) * MDP(25) + (-t104 * t128 + t114 * t123) * MDP(26) + (t100 * t121 + t101 * t119 + t104 * t114) * MDP(27) + t207 * t187 + (MDP(7) + (t150 - t152) * MDP(16)) * t159 + (-MDP(10) + MDP(13)) * t148 + ((-MDP(9) + MDP(12)) * pkin(7) - t207 * pkin(2) + t158 * MDP(17) - t156 * MDP(18) + MDP(6) + t204) * t157; MDP(8) + t152 * MDP(15) - 0.2e1 * t178 + (t145 ^ 2 + t180) * MDP(23) + t180 * MDP(27) + (t183 - 0.2e1 * t184 + 0.2e1 * t185) * t114 + (-0.2e1 * MDP(12) + t193) * pkin(2) + (MDP(14) * qJ(3) + t156 * t201 + t158 * t200 + 0.2e1 * MDP(13)) * qJ(3) - (t199 + t197) * t208; t168 * MDP(23) + t169 * MDP(27) + (t162 + t207) * t157 + t206 * (t122 * t165 + t123 * t128); t205 * t208 - t206 * t209 + MDP(12) - t193; t205 * t209 + MDP(14); t116 * MDP(20) - t117 * MDP(21) + t102 * MDP(24) + (t122 * t139 + t123 * t143) * MDP(25) + t100 * MDP(26) + (t100 * t139 - t101 * t143) * MDP(27) + t164 * t159 + (MDP(19) + (pkin(5) + t143) * MDP(24) + t182) * t157 + ((t122 * t154 + t123 * t155) * MDP(22) + (t102 * t155 + t103 * t154) * MDP(23)) * pkin(4); -t167 * MDP(25) + (-t119 * t143 + t121 * t139) * MDP(27) + (-MDP(20) * t196 + MDP(17)) * t158 + (MDP(21) * t196 - MDP(18)) * t156 + (-t166 * MDP(22) + (-t119 * t155 + t121 * t154) * MDP(23)) * pkin(4) + t204; MDP(23) * pkin(4) * t166 + MDP(27) * t167 + t162; MDP(19) + (t139 ^ 2 + t143 ^ 2) * MDP(27) + (t154 ^ 2 + t155 ^ 2) * MDP(23) * pkin(4) ^ 2 + t143 * t198 + 0.2e1 * t182; MDP(23) * t124 - MDP(24) * t122 + MDP(26) * t123 + MDP(27) * t104; MDP(23) * t145 + t183 - t184 + t185; 0; 0; t205; -MDP(24) * t157 - MDP(25) * t123 + MDP(27) * t101; MDP(25) * t128 + MDP(27) * t119; -t128 * MDP(27); -MDP(27) * t143 - MDP(24); 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
