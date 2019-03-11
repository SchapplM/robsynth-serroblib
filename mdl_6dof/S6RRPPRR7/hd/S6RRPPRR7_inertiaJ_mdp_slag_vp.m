% Calculate joint inertia matrix for
% S6RRPPRR7
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
%   see S6RRPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR7_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:21:01
% EndTime: 2019-03-09 09:21:04
% DurationCPUTime: 0.94s
% Computational Cost: add. (679->196), mult. (1427->268), div. (0->0), fcn. (1382->8), ass. (0->84)
t145 = sin(pkin(6));
t153 = cos(qJ(2));
t192 = t145 * t153;
t150 = sin(qJ(2));
t193 = t145 * t150;
t117 = -t145 * pkin(1) - pkin(2) * t192 - qJ(3) * t193;
t113 = pkin(3) * t192 - t117;
t106 = (pkin(4) * t150 + pkin(9) * t153) * t145 + t113;
t154 = -pkin(2) - pkin(3);
t139 = -pkin(9) + t154;
t146 = cos(pkin(6));
t198 = pkin(1) * t146;
t188 = pkin(8) * t193 - t153 * t198;
t160 = -qJ(4) * t193 + t188;
t107 = t139 * t146 + t160;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t103 = t152 * t106 - t149 * t107;
t104 = t149 * t106 + t152 * t107;
t121 = -t146 * t149 - t152 * t192;
t207 = t121 * MDP(21) + t103 * MDP(24) - t104 * MDP(25);
t148 = sin(qJ(6));
t151 = cos(qJ(6));
t161 = t148 * MDP(31) + t151 * MDP(32);
t157 = t148 * MDP(28) + t151 * MDP(29) - pkin(10) * t161;
t206 = MDP(22) - t157;
t205 = 0.2e1 * t146;
t173 = -MDP(11) + MDP(16);
t203 = MDP(15) + MDP(13);
t202 = 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * t154 * MDP(16) + MDP(8);
t200 = 0.2e1 * MDP(31);
t199 = 0.2e1 * MDP(32);
t147 = qJ(3) + pkin(4);
t101 = -pkin(5) * t193 - t103;
t197 = t101 * t148;
t196 = t101 * t151;
t195 = t139 * t148;
t194 = t139 * t151;
t191 = t148 * t152;
t119 = -t146 * t152 + t149 * t192;
t190 = t149 * t119;
t189 = t151 * t152;
t123 = pkin(8) * t192 + t150 * t198;
t110 = t121 * t148 - t151 * t193;
t186 = t110 * MDP(29);
t111 = t121 * t151 + t148 * t193;
t185 = t111 * MDP(26);
t184 = t111 * MDP(28);
t183 = t119 * MDP(30);
t182 = t121 * MDP(20);
t180 = t147 * MDP(25);
t178 = t149 * MDP(25);
t177 = t151 * MDP(26);
t175 = t152 * MDP(24);
t174 = t152 * MDP(25);
t172 = MDP(12) - MDP(17);
t135 = t146 * qJ(3);
t116 = t135 + t123;
t170 = t148 * t151 * MDP(27);
t168 = -t123 * MDP(10) - t188 * MDP(9);
t166 = -t149 * MDP(24) - t174;
t165 = -MDP(28) * t151 + MDP(29) * t148;
t124 = t152 * pkin(5) + t149 * pkin(10) + t147;
t114 = t151 * t124 - t139 * t191;
t115 = t148 * t124 + t139 * t189;
t163 = t114 * MDP(31) - t115 * MDP(32);
t162 = t151 * MDP(31) - t148 * MDP(32);
t159 = -MDP(20) - t165;
t158 = MDP(24) + t162;
t112 = qJ(4) * t192 - t116;
t108 = t146 * pkin(4) - t112;
t102 = pkin(10) * t193 + t104;
t105 = -t119 * pkin(5) - t121 * pkin(10) + t108;
t100 = t151 * t102 + t148 * t105;
t99 = -t148 * t102 + t151 * t105;
t156 = t99 * MDP(31) - t100 * MDP(32) - t183 + t184 - t186;
t155 = qJ(3) ^ 2;
t144 = t152 ^ 2;
t143 = t151 ^ 2;
t141 = t149 ^ 2;
t140 = t148 ^ 2;
t118 = -t146 * pkin(2) + t188;
t109 = t154 * t146 + t160;
t1 = [t146 ^ 2 * MDP(8) + (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) * MDP(14) + t121 ^ 2 * MDP(19) + (t109 ^ 2 + t112 ^ 2 + t113 ^ 2) * MDP(18) + MDP(1) + (-0.2e1 * t110 * MDP(27) + t185) * t111 + (0.2e1 * MDP(22) * t193 + 0.2e1 * t182 + t183 - 0.2e1 * t184 + 0.2e1 * t186) * t119 + (t100 * t119 + t101 * t111) * t199 + (t101 * t110 - t99 * t119) * t200 + 0.2e1 * (-t119 * MDP(24) + t121 * MDP(25)) * t108 + (-t118 * MDP(11) + t116 * MDP(13) - t112 * MDP(15) + t109 * MDP(16) + t168) * t205 + ((t150 * MDP(6) + t153 * MDP(7)) * t205 + 0.2e1 * (-t117 * MDP(11) + t116 * MDP(12) - t113 * MDP(16) + t112 * MDP(17)) * t153 + ((MDP(23) + MDP(4)) * t150 ^ 2 + 0.2e1 * (-MDP(10) * t150 + MDP(9) * t153) * pkin(1)) * t145 + 0.2e1 * (t118 * MDP(12) - t117 * MDP(13) + t113 * MDP(15) - t109 * MDP(17) + MDP(5) * t192 + t207) * t150) * t145; (-t112 * qJ(3) + t109 * t154) * MDP(18) + t121 * t180 + (-t118 * pkin(2) + t116 * qJ(3)) * MDP(14) + t202 * t146 + (-t147 * MDP(24) - t163) * t119 + (t108 * MDP(24) + t156 - t182) * t152 + ((t110 * t151 + t111 * t148) * MDP(27) + (t111 * t139 - t196) * MDP(32) + (t110 * t139 - t197) * MDP(31) - t108 * MDP(25) - t111 * t177 - t121 * MDP(19) + t159 * t119) * t149 + ((-qJ(4) * MDP(15) + t172 * qJ(3) + MDP(7)) * t153 + (-pkin(2) * MDP(12) - qJ(4) * MDP(16) - t154 * MDP(17) - t149 * MDP(21) - t152 * MDP(22) + t166 * t139 + MDP(6)) * t150) * t145 + t168 + t173 * t188 + t203 * (0.2e1 * t135 + t123); (pkin(2) ^ 2 + t155) * MDP(14) + (t154 ^ 2 + t155) * MDP(18) + 0.2e1 * t147 * t175 + t144 * MDP(30) + (t143 * MDP(26) + MDP(19) - 0.2e1 * t170) * t141 + (t114 * t152 - t141 * t195) * t200 + (-t115 * t152 - t141 * t194) * t199 + 0.2e1 * t203 * qJ(3) + t202 - 0.2e1 * (t152 * t159 + t180) * t149; t118 * MDP(14) + t109 * MDP(18) + (t149 * t110 + t119 * t191) * MDP(31) + (t149 * t111 + t119 * t189) * MDP(32) + t173 * t146 + (t166 + t172) * t193; -pkin(2) * MDP(14) + t154 * MDP(18) + t173 + t161 * (-t141 - t144); MDP(14) + MDP(18); t113 * MDP(18) + (-t152 * t110 + t148 * t190) * MDP(31) + (-t111 * t152 + t151 * t190) * MDP(32) + (-t153 * MDP(16) + (MDP(15) + t175 - t178) * t150) * t145; 0; 0; MDP(18); MDP(23) * t193 + t148 * t185 + (-t148 * t110 + t111 * t151) * MDP(27) + (-pkin(5) * t110 - t196) * MDP(31) + (-pkin(5) * t111 + t197) * MDP(32) + t206 * t119 + t207; (-t139 * MDP(25) - t206) * t152 + (-MDP(21) - t139 * MDP(24) - t148 * t177 + (t140 - t143) * MDP(27) + (pkin(5) * t148 - t194) * MDP(31) + (pkin(5) * t151 + t195) * MDP(32)) * t149; -t149 * t158 - t174; t152 * t158 - t178; t140 * MDP(26) + 0.2e1 * pkin(5) * t162 + MDP(23) + 0.2e1 * t170; t156; t152 * MDP(30) + t149 * t165 + t163; -t161 * t152; -t161 * t149; t157; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
