% Calculate joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR8_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:22
% EndTime: 2019-03-09 10:53:25
% DurationCPUTime: 0.83s
% Computational Cost: add. (1023->194), mult. (1999->277), div. (0->0), fcn. (2080->8), ass. (0->83)
t150 = sin(qJ(6));
t153 = cos(qJ(6));
t152 = sin(qJ(2));
t148 = sin(pkin(10));
t149 = cos(pkin(10));
t151 = sin(qJ(4));
t154 = cos(qJ(4));
t202 = -t151 * t148 + t154 * t149;
t123 = t202 * t152;
t155 = cos(qJ(2));
t143 = t155 * pkin(4);
t133 = -t155 * pkin(2) - t152 * qJ(3) - pkin(1);
t126 = t149 * t133;
t192 = pkin(7) * t148;
t113 = -t149 * t152 * pkin(8) + t126 + (-pkin(3) - t192) * t155;
t191 = pkin(7) * t155;
t120 = t148 * t133 + t149 * t191;
t187 = t148 * t152;
t116 = -pkin(8) * t187 + t120;
t183 = -t154 * t113 + t151 * t116;
t98 = t143 + t183;
t93 = t155 * pkin(5) - t123 * pkin(9) + t98;
t129 = t154 * t148 + t151 * t149;
t122 = t129 * t152;
t101 = t151 * t113 + t154 * t116;
t97 = -t155 * qJ(5) + t101;
t94 = t122 * pkin(9) + t97;
t169 = t150 * t94 - t153 * t93;
t92 = t150 * t93 + t153 * t94;
t203 = t169 * MDP(31) + t92 * MDP(32);
t201 = MDP(14) * qJ(3);
t171 = MDP(21) - MDP(24);
t200 = 2 * MDP(13);
t199 = -2 * MDP(16);
t198 = 2 * MDP(21);
t197 = 2 * MDP(22);
t196 = 2 * MDP(23);
t195 = -2 * MDP(27);
t194 = 0.2e1 * MDP(32);
t193 = -pkin(4) - pkin(5);
t190 = pkin(8) + qJ(3);
t189 = pkin(2) * MDP(14);
t188 = pkin(4) * MDP(25);
t130 = pkin(3) * t187 + t152 * pkin(7);
t181 = qJ(5) * MDP(24);
t106 = t150 * t122 + t153 * t123;
t180 = t106 * MDP(26);
t111 = t150 * t129 + t153 * t202;
t179 = t111 * MDP(31);
t178 = t123 * MDP(15);
t177 = (t150 * qJ(5) - t153 * t193) * MDP(31);
t176 = (t153 * qJ(5) + t150 * t193) * MDP(32);
t175 = t148 * MDP(12);
t174 = t149 * MDP(11);
t173 = MDP(19) + MDP(30);
t172 = MDP(20) + MDP(22);
t140 = -t149 * pkin(3) - pkin(2);
t134 = t190 * t148;
t135 = t190 * t149;
t117 = t154 * t134 + t151 * t135;
t170 = t123 * qJ(5) - t130;
t118 = -t151 * t134 + t154 * t135;
t167 = t148 * MDP(11) + t149 * MDP(12);
t105 = -t153 * t122 + t150 * t123;
t166 = t106 * MDP(28) - t105 * MDP(29);
t165 = t153 * MDP(31) - t150 * MDP(32);
t164 = t129 * qJ(5) - t140;
t163 = -t129 * pkin(9) + t117;
t162 = MDP(22) + t165;
t161 = -MDP(30) - t176 - t177;
t160 = -t174 + t175 - t189;
t107 = -pkin(9) * t202 + t118;
t112 = t153 * t129 - t150 * t202;
t159 = t112 * MDP(28) - t111 * MDP(29) - (t150 * t107 - t153 * t163) * MDP(31) - (t153 * t107 + t150 * t163) * MDP(32);
t158 = t123 * MDP(17) - t122 * MDP(18) - t166;
t157 = t129 * MDP(17) + MDP(18) * t202 - t159;
t146 = t152 ^ 2;
t119 = -t148 * t191 + t126;
t110 = -pkin(4) * t202 - t164;
t103 = -t193 * t202 + t164;
t102 = t122 * pkin(4) - t170;
t99 = t193 * t122 + t170;
t1 = [(t102 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(25) + (t146 * pkin(7) ^ 2 + t119 ^ 2 + t120 ^ 2) * MDP(14) + t146 * MDP(4) - 0.2e1 * pkin(1) * t152 * MDP(10) + MDP(1) + t173 * t155 ^ 2 + (t122 * t199 + t178) * t123 + (t105 * t195 + t180) * t106 + 0.2e1 * (t152 * MDP(5) + pkin(1) * MDP(9) - t158) * t155 + 0.2e1 * (-t119 * t155 + t146 * t192) * MDP(11) + 0.2e1 * (t146 * pkin(7) * t149 + t120 * t155) * MDP(12) + 0.2e1 * (t130 * t122 + t155 * t183) * MDP(20) + (t101 * t155 + t130 * t123) * t198 + (t102 * t122 + t98 * t155) * t197 + 0.2e1 * (-t102 * t123 - t97 * t155) * MDP(24) + (t99 * t106 - t92 * t155) * t194 + 0.2e1 * (t99 * t105 - t155 * t169) * MDP(31) + (-t97 * t122 + t98 * t123) * t196 + (-t119 * t149 - t120 * t148) * t152 * t200; (-t129 * t122 + t123 * t202) * MDP(16) + (t140 * t122 - t130 * t202) * MDP(20) + (t140 * t123 + t130 * t129) * MDP(21) + (-t102 * t202 + t110 * t122) * MDP(22) + (t117 * t123 - t118 * t122 + t98 * t129 + t202 * t97) * MDP(23) + (-t102 * t129 - t110 * t123) * MDP(24) + (t102 * t110 + t98 * t117 + t97 * t118) * MDP(25) + (-t112 * t105 - t106 * t111) * MDP(27) + (t103 * t105 + t99 * t111) * MDP(31) + (t103 * t106 + t99 * t112) * MDP(32) + t129 * t178 + t112 * t180 + (MDP(6) - t167 * pkin(2) + (-MDP(9) + t160) * pkin(7)) * t152 + (-pkin(7) * MDP(10) + t167 * qJ(3) + t172 * t117 + t171 * t118 + MDP(7) - t157) * t155 + (MDP(13) + t201) * (-t119 * t148 + t120 * t149); MDP(8) + (t110 ^ 2 + t117 ^ 2 + t118 ^ 2) * MDP(25) + 0.2e1 * t103 * t179 + (0.2e1 * t174 - 0.2e1 * t175 + t189) * pkin(2) - 0.2e1 * (t140 * MDP(20) + t110 * MDP(22) - t118 * MDP(23)) * t202 + (MDP(15) * t129 - 0.2e1 * t110 * MDP(24) + t117 * t196 + t140 * t198 - t199 * t202) * t129 + (MDP(26) * t112 + t103 * t194 + t111 * t195) * t112 + (t200 + t201) * (t148 ^ 2 + t149 ^ 2) * qJ(3); t102 * MDP(25) - t105 * MDP(31) - t106 * MDP(32) + t171 * t123 + t172 * t122 + (pkin(7) * MDP(14) + t167) * t152; t110 * MDP(25) - t112 * MDP(32) + t171 * t129 - t172 * t202 + t160 - t179; MDP(14) + MDP(25); -t183 * MDP(20) + (-0.2e1 * t143 - t183) * MDP(22) + (-pkin(4) * t123 - t122 * qJ(5)) * MDP(23) + (-t98 * pkin(4) + t97 * qJ(5)) * MDP(25) + (-MDP(19) + t161 - 0.2e1 * t181) * t155 + t158 - t171 * t101 + t203; (-pkin(4) * t129 + qJ(5) * t202) * MDP(23) + (MDP(25) * qJ(5) - t171) * t118 + (-t172 - t188) * t117 + t157; 0; pkin(4) * t197 + 0.2e1 * t181 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(25) + 0.2e1 * t177 + 0.2e1 * t176 + t173; t123 * MDP(23) + t98 * MDP(25) + t162 * t155; t129 * MDP(23) + t117 * MDP(25); 0; -t162 - t188; MDP(25); t155 * MDP(30) + t166 - t203; t159; 0; t161; t165; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
