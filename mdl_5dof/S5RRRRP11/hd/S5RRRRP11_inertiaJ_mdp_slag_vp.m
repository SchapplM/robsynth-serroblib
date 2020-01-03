% Calculate joint inertia matrix for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP11_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:35
% EndTime: 2019-12-31 22:18:38
% DurationCPUTime: 0.94s
% Computational Cost: add. (1024->222), mult. (2292->319), div. (0->0), fcn. (2337->8), ass. (0->90)
t139 = sin(pkin(5));
t198 = 0.2e1 * t139;
t142 = sin(qJ(3));
t197 = -0.2e1 * t142;
t196 = 2 * MDP(25);
t195 = 2 * MDP(26);
t143 = sin(qJ(2));
t194 = pkin(1) * t143;
t146 = cos(qJ(2));
t193 = pkin(1) * t146;
t141 = sin(qJ(4));
t192 = pkin(8) * t141;
t144 = cos(qJ(4));
t191 = pkin(8) * t144;
t140 = cos(pkin(5));
t145 = cos(qJ(3));
t182 = t139 * t143;
t121 = -t140 * t145 + t142 * t182;
t190 = t121 * pkin(4);
t189 = pkin(2) * MDP(16);
t188 = pkin(2) * MDP(17);
t181 = t139 * t146;
t164 = pkin(7) * t181;
t118 = t164 + (pkin(8) + t194) * t140;
t119 = (-pkin(2) * t146 - pkin(8) * t143 - pkin(1)) * t139;
t108 = -t142 * t118 + t145 * t119;
t106 = pkin(3) * t181 - t108;
t187 = t106 * t141;
t186 = t106 * t144;
t122 = t140 * t142 + t145 * t182;
t110 = t122 * t141 + t144 * t181;
t185 = t110 * t144;
t111 = t122 * t144 - t141 * t181;
t184 = t111 * t141;
t183 = t121 * qJ(5);
t180 = t140 * MDP(8);
t179 = t141 * t145;
t178 = t143 * MDP(6);
t130 = pkin(7) * t182;
t117 = t130 + (-pkin(2) - t193) * t140;
t105 = t121 * pkin(3) - t122 * pkin(9) + t117;
t109 = t145 * t118 + t142 * t119;
t107 = -pkin(9) * t181 + t109;
t100 = t141 * t105 + t144 * t107;
t128 = -t145 * pkin(3) - t142 * pkin(9) - pkin(2);
t116 = t141 * t128 + t145 * t191;
t136 = t141 ^ 2;
t138 = t144 ^ 2;
t177 = t136 + t138;
t176 = MDP(15) * t146;
t101 = t110 * pkin(4) - t111 * qJ(5) + t106;
t175 = t101 * MDP(28);
t174 = t110 * MDP(21);
t173 = t111 * MDP(18);
t172 = t111 * MDP(20);
t171 = t121 * MDP(22);
t170 = t122 * MDP(12);
t169 = t122 * MDP(13);
t168 = t141 * MDP(20);
t167 = t144 * MDP(18);
t166 = t144 * MDP(19);
t165 = MDP(24) - MDP(27);
t163 = pkin(8) * MDP(16) - MDP(13);
t162 = pkin(8) * MDP(17) - MDP(14);
t161 = -MDP(28) * pkin(4) - MDP(25);
t160 = -t144 * t105 + t141 * t107;
t159 = -0.2e1 * qJ(5) * MDP(27) - MDP(22);
t158 = -t144 * pkin(4) - t141 * qJ(5);
t157 = -pkin(4) * t141 + t144 * qJ(5);
t97 = t100 + t183;
t98 = t160 - t190;
t156 = t98 * t141 + t97 * t144;
t112 = -t145 * qJ(5) + t116;
t126 = t144 * t128;
t113 = -t126 + (pkin(4) + t192) * t145;
t155 = t112 * t144 + t113 * t141;
t154 = t144 * MDP(20) - t141 * MDP(21);
t153 = t144 * MDP(21) + t168;
t152 = (-pkin(8) * t179 + t126) * MDP(23) - t116 * MDP(24);
t151 = t141 * MDP(25) - t144 * MDP(27);
t150 = -MDP(12) + t154;
t149 = t113 * MDP(25) - t112 * MDP(27) - t152;
t148 = -MDP(23) * t160 - t100 * MDP(24) + t171 + t172 - t174;
t135 = t139 ^ 2;
t132 = pkin(9) * t179;
t127 = -pkin(3) + t158;
t124 = t140 * t194 + t164;
t123 = t140 * t193 - t130;
t120 = (pkin(8) - t157) * t142;
t1 = [t122 ^ 2 * MDP(11) + (t101 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(28) + t135 * t143 ^ 2 * MDP(4) + MDP(1) + (t178 * t198 + t180) * t140 + (-0.2e1 * t110 * MDP(19) + t173) * t111 + ((MDP(7) * t140 - t169) * t198 + (0.2e1 * MDP(5) * t143 + t176) * t135) * t146 + (0.2e1 * MDP(14) * t181 - 0.2e1 * t170 + t171 + 0.2e1 * t172 - 0.2e1 * t174) * t121 + (-t97 * t110 + t98 * t111) * t195 + 0.2e1 * (t106 * t110 - t121 * t160) * MDP(23) + 0.2e1 * (-t101 * t111 + t97 * t121) * MDP(27) + (t101 * t110 - t98 * t121) * t196 + 0.2e1 * (-t100 * t121 + t106 * t111) * MDP(24) + 0.2e1 * (-t124 * t140 - t135 * t194) * MDP(10) + 0.2e1 * (-t108 * t181 + t117 * t121) * MDP(16) + 0.2e1 * (t109 * t181 + t117 * t122) * MDP(17) + 0.2e1 * (t123 * t140 + t135 * t193) * MDP(9); t180 + t123 * MDP(9) - t124 * MDP(10) - t122 * t188 + (-t112 * t110 + t113 * t111) * MDP(26) + (t97 * t112 + t98 * t113) * MDP(28) + (t146 * MDP(7) + t178) * t139 + (t110 * MDP(25) - t111 * MDP(27) + t175) * t120 + (-t149 - t189) * t121 + (-t117 * MDP(16) + t98 * MDP(25) - t97 * MDP(27) + t162 * t181 - t148 + t170) * t145 + (t122 * MDP(11) + t117 * MDP(17) + t111 * t167 + (-t184 - t185) * MDP(19) + (pkin(8) * t110 + t187) * MDP(23) + (pkin(8) * t111 + t186) * MDP(24) + (-t141 * t97 + t144 * t98) * MDP(26) + t163 * t181 + t151 * t101 + t150 * t121) * t142; MDP(8) + t188 * t197 + (t112 ^ 2 + t113 ^ 2 + t120 ^ 2) * MDP(28) + 0.2e1 * ((-t112 * t141 + t113 * t144) * MDP(26) + t151 * t120) * t142 + (t145 * MDP(22) + t150 * t197 + 0.2e1 * t149 + 0.2e1 * t189) * t145 + (t138 * MDP(18) - 0.2e1 * t141 * t166 + MDP(11) + 0.2e1 * (t141 * MDP(23) + t144 * MDP(24)) * pkin(8)) * t142 ^ 2; t169 - t139 * t176 + t108 * MDP(16) - t109 * MDP(17) + t141 * t173 + (-t141 * t110 + t111 * t144) * MDP(19) + (-pkin(3) * t110 - t186) * MDP(23) + (-pkin(3) * t111 + t187) * MDP(24) + (-t101 * t144 + t127 * t110) * MDP(25) + t156 * MDP(26) + (-t101 * t141 - t127 * t111) * MDP(27) + t127 * t175 + ((t184 - t185) * MDP(26) + t156 * MDP(28)) * pkin(9) + (-MDP(14) + (-t165 * t144 + (-MDP(23) - MDP(25)) * t141) * pkin(9) + t153) * t121; t132 * MDP(23) + (-t120 * t144 + t132) * MDP(25) + t155 * MDP(26) - t120 * t141 * MDP(27) + (t155 * pkin(9) + t120 * t127) * MDP(28) + (-t168 + (t165 * pkin(9) - MDP(21)) * t144 - t162) * t145 + (t141 * t167 + (-t136 + t138) * MDP(19) + (-pkin(3) * t141 - t191) * MDP(23) + (-pkin(3) * t144 + t192) * MDP(24) + t151 * t127 - t163) * t142; MDP(15) + t136 * MDP(18) + (t177 * pkin(9) ^ 2 + t127 ^ 2) * MDP(28) + t177 * pkin(9) * t195 + 0.2e1 * (pkin(3) * MDP(23) - t127 * MDP(25)) * t144 + 0.2e1 * (-pkin(3) * MDP(24) - t127 * MDP(27) + t166) * t141; (-t160 + 0.2e1 * t190) * MDP(25) + (-t111 * pkin(4) - t110 * qJ(5)) * MDP(26) + (t100 + 0.2e1 * t183) * MDP(27) + (-t98 * pkin(4) + t97 * qJ(5)) * MDP(28) + t148; t126 * MDP(25) + t116 * MDP(27) + (-t113 * pkin(4) + t112 * qJ(5)) * MDP(28) + ((-0.2e1 * pkin(4) - t192) * MDP(25) + t159) * t145 + (t158 * MDP(26) + t154) * t142 + t152; t157 * MDP(26) + ((MDP(28) * qJ(5) - t165) * t144 + (-MDP(23) + t161) * t141) * pkin(9) + t153; pkin(4) * t196 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(28) - t159; -t121 * MDP(25) + t111 * MDP(26) + t98 * MDP(28); t144 * t142 * MDP(26) + t145 * MDP(25) + t113 * MDP(28); (MDP(28) * pkin(9) + MDP(26)) * t141; t161; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
