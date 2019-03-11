% Calculate joint inertia matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP11_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:38
% EndTime: 2019-03-09 12:48:40
% DurationCPUTime: 0.56s
% Computational Cost: add. (721->162), mult. (1240->225), div. (0->0), fcn. (1199->6), ass. (0->77)
t182 = pkin(3) + pkin(7);
t137 = sin(qJ(4));
t140 = cos(qJ(4));
t142 = -pkin(2) - pkin(8);
t172 = -pkin(9) + t142;
t118 = t172 * t137;
t119 = t172 * t140;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t100 = t118 * t139 + t119 * t136;
t115 = t136 * t140 + t137 * t139;
t116 = -t136 * t137 + t139 * t140;
t99 = -t118 * t136 + t139 * t119;
t151 = t116 * MDP(24) - t115 * MDP(25) + t99 * MDP(27) - t100 * MDP(28);
t181 = (MDP(20) * t142 + MDP(17)) * t140 + (-MDP(21) * t142 - MDP(18)) * t137 + t151;
t141 = cos(qJ(2));
t165 = t140 * t141;
t166 = t137 * t141;
t104 = t136 * t166 - t139 * t165;
t105 = t115 * t141;
t180 = -t105 * MDP(24) + t104 * MDP(25);
t114 = t116 ^ 2;
t179 = t115 ^ 2 + t114;
t177 = 0.2e1 * t141;
t176 = -2 * MDP(23);
t175 = pkin(4) * t136;
t138 = sin(qJ(2));
t174 = pkin(4) * t138;
t173 = pkin(4) * t139;
t171 = MDP(30) * pkin(5);
t170 = pkin(2) * MDP(14);
t154 = -qJ(3) * t138 - pkin(1);
t112 = t141 * t142 + t154;
t153 = pkin(9) * t141 - t112;
t121 = t182 * t138;
t168 = t121 * t137;
t93 = -t140 * t153 + t168;
t169 = t139 * t93;
t167 = t137 * t140;
t125 = t137 * pkin(4) + qJ(3);
t164 = t116 * MDP(27) - t115 * MDP(28);
t122 = t182 * t141;
t163 = t105 * MDP(22);
t161 = t137 * MDP(20);
t160 = t139 * MDP(27);
t159 = t140 * MDP(21);
t158 = pkin(7) ^ 2 * MDP(14);
t157 = MDP(19) + MDP(26);
t156 = t138 * MDP(26) + t180;
t107 = pkin(4) * t165 + t122;
t155 = MDP(16) * t167;
t117 = t140 * t121;
t92 = t137 * t153 + t117 + t174;
t87 = -t136 * t93 + t139 * t92;
t152 = MDP(12) - t170;
t85 = pkin(5) * t138 + qJ(6) * t105 + t87;
t88 = t136 * t92 + t169;
t86 = qJ(6) * t104 + t88;
t150 = t115 * t86 + t116 * t85;
t90 = -qJ(6) * t116 + t99;
t91 = -qJ(6) * t115 + t100;
t149 = t115 * t91 + t116 * t90;
t148 = -MDP(17) * t137 - MDP(18) * t140;
t147 = t140 * MDP(20) - t137 * MDP(21);
t127 = pkin(5) + t173;
t146 = t115 * t175 + t116 * t127;
t145 = pkin(7) * MDP(14) + MDP(11) + t147;
t135 = t141 ^ 2;
t134 = t140 ^ 2;
t133 = t138 ^ 2;
t132 = t137 ^ 2;
t120 = -pkin(2) * t141 + t154;
t101 = pkin(5) * t115 + t125;
t96 = t112 * t140 + t168;
t95 = -t112 * t137 + t117;
t94 = -pkin(5) * t104 + t107;
t1 = [(t85 ^ 2 + t86 ^ 2 + t94 ^ 2) * MDP(30) + pkin(1) * MDP(9) * t177 + MDP(1) + (MDP(12) * t177 + MDP(14) * t120) * t120 + (t132 * MDP(15) + 0.2e1 * t155 + t158) * t135 + (t104 * t176 + t163) * t105 + (MDP(4) + t157 + t158) * t133 + 0.2e1 * (-pkin(1) * MDP(10) - t120 * MDP(13) + (MDP(5) + t148) * t141 + t180) * t138 + 0.2e1 * (t104 * t86 + t105 * t85) * MDP(29) + 0.2e1 * (-t104 * t107 + t138 * t87) * MDP(27) + 0.2e1 * (-t105 * t107 - t138 * t88) * MDP(28) + 0.2e1 * (-t122 * t166 - t138 * t96) * MDP(21) + 0.2e1 * (t122 * t165 + t138 * t95) * MDP(20) + 0.2e1 * (t133 + t135) * MDP(11) * pkin(7); -t116 * t163 + (t104 * t116 + t105 * t115) * MDP(23) + (-t104 * t125 + t107 * t115) * MDP(27) + (-t105 * t125 + t107 * t116) * MDP(28) + (t104 * t91 + t105 * t90 - t150) * MDP(29) + (t101 * t94 + t85 * t90 + t86 * t91) * MDP(30) + (t159 + t161) * t122 + (MDP(7) - MDP(15) * t167 + (t132 - t134) * MDP(16) + (-MDP(10) + MDP(13)) * pkin(7) + t145 * qJ(3)) * t141 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(9) + t152) * pkin(7) + t181) * t138; MDP(8) + t134 * MDP(15) - 0.2e1 * t155 + t114 * MDP(22) + t116 * t115 * t176 - 0.2e1 * t149 * MDP(29) + (t101 ^ 2 + t90 ^ 2 + t91 ^ 2) * MDP(30) + 0.2e1 * (t115 * MDP(27) + t116 * MDP(28)) * t125 + (-0.2e1 * MDP(12) + t170) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * t159 + 0.2e1 * t161) * qJ(3); (t104 * t115 + t105 * t116) * MDP(29) + t150 * MDP(30) + (t145 + t164) * t138; -t179 * MDP(29) + t149 * MDP(30) + t152; t179 * MDP(30) + MDP(14); t138 * MDP(19) + t95 * MDP(20) - t96 * MDP(21) + (t138 * t173 + t87) * MDP(27) + (-t169 + (-t92 - t174) * t136) * MDP(28) + (t104 * t175 + t105 * t127) * MDP(29) + (t127 * t85 + t175 * t86) * MDP(30) + t148 * t141 + t156; -t146 * MDP(29) + (t127 * t90 + t175 * t91) * MDP(30) + t181; MDP(30) * t146 + t147 + t164; t127 ^ 2 * MDP(30) + (0.2e1 * t160 + (MDP(30) * t175 - 0.2e1 * MDP(28)) * t136) * pkin(4) + t157; t87 * MDP(27) - t88 * MDP(28) + (t105 * MDP(29) + t85 * MDP(30)) * pkin(5) + t156; (-t116 * MDP(29) + MDP(30) * t90) * pkin(5) + t151; t116 * t171 + t164; t127 * t171 + MDP(26) + (-MDP(28) * t136 + t160) * pkin(4); MDP(30) * pkin(5) ^ 2 + MDP(26); t94 * MDP(30); t101 * MDP(30); 0; 0; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
