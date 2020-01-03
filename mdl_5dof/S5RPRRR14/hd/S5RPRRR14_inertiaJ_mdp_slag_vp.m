% Calculate joint inertia matrix for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR14_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:46
% EndTime: 2019-12-31 19:19:48
% DurationCPUTime: 0.69s
% Computational Cost: add. (1216->173), mult. (3303->275), div. (0->0), fcn. (3749->12), ass. (0->88)
t144 = sin(pkin(5));
t143 = sin(pkin(6));
t183 = cos(pkin(5));
t164 = t183 * t143;
t145 = cos(pkin(11));
t146 = cos(pkin(6));
t177 = t145 * t146;
t194 = t144 * t177 + t164;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t154 = -(MDP(27) * t147 + MDP(28) * t150) * pkin(10) + t147 * MDP(24) + t150 * MDP(25);
t142 = sin(pkin(11));
t166 = pkin(1) * t183;
t178 = t144 * t145;
t127 = qJ(2) * t178 + t142 * t166;
t113 = t194 * pkin(8) + t127;
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t135 = t145 * t166;
t181 = t142 * t144;
t116 = t183 * pkin(2) + t135 + (-pkin(8) * t146 - qJ(2)) * t181;
t120 = (-pkin(8) * t142 * t143 - pkin(2) * t145 - pkin(1)) * t144;
t161 = t116 * t146 + t120 * t143;
t104 = -t149 * t113 + t161 * t152;
t193 = 0.2e1 * MDP(27);
t192 = 0.2e1 * MDP(28);
t138 = t144 ^ 2;
t191 = pkin(1) * t138;
t190 = pkin(9) * t147;
t189 = pkin(9) * t150;
t151 = cos(qJ(4));
t188 = pkin(9) * t151;
t187 = MDP(20) * pkin(3);
t186 = pkin(3) * MDP(21);
t114 = t149 * t181 - t194 * t152;
t108 = -t116 * t143 + t146 * t120;
t115 = t149 * t164 + (t142 * t152 + t149 * t177) * t144;
t101 = pkin(3) * t114 - pkin(9) * t115 + t108;
t105 = t152 * t113 + t161 * t149;
t124 = t143 * t178 - t183 * t146;
t103 = -t124 * pkin(9) + t105;
t148 = sin(qJ(4));
t98 = t101 * t151 - t103 * t148;
t96 = -pkin(4) * t114 - t98;
t185 = t147 * t96;
t184 = t150 * t96;
t180 = t143 * t149;
t128 = -t151 * t146 + t148 * t180;
t182 = t128 * t148;
t179 = t143 * t152;
t110 = t115 * t151 - t124 * t148;
t107 = t110 * t150 + t114 * t147;
t176 = MDP(22) * t107;
t175 = MDP(22) * t150;
t106 = t110 * t147 - t114 * t150;
t174 = t106 * MDP(25);
t173 = t107 * MDP(24);
t109 = t115 * t148 + t124 * t151;
t172 = t109 * MDP(26);
t171 = t110 * MDP(16);
t170 = t110 * MDP(17);
t169 = t114 * MDP(19);
t168 = t151 * MDP(26);
t165 = t147 * t150 * MDP(23);
t163 = -pkin(9) * MDP(20) + MDP(17);
t162 = -pkin(9) * MDP(21) + MDP(18);
t99 = t101 * t148 + t103 * t151;
t160 = MDP(24) * t150 - MDP(25) * t147;
t131 = -t151 * pkin(4) - t148 * pkin(10) - pkin(3);
t121 = t150 * t131 - t147 * t188;
t122 = t147 * t131 + t150 * t188;
t158 = t121 * MDP(27) - t122 * MDP(28);
t157 = MDP(27) * t150 - MDP(28) * t147;
t155 = -MDP(16) + t160;
t102 = t124 * pkin(3) - t104;
t100 = t109 * pkin(4) - t110 * pkin(10) + t102;
t97 = pkin(10) * t114 + t99;
t94 = t100 * t150 - t147 * t97;
t95 = t100 * t147 + t150 * t97;
t153 = t94 * MDP(27) - t95 * MDP(28) + t172 + t173 - t174;
t141 = t150 ^ 2;
t140 = t148 ^ 2;
t139 = t147 ^ 2;
t129 = t148 * t146 + t151 * t180;
t126 = -qJ(2) * t181 + t135;
t118 = t150 * t129 - t147 * t179;
t117 = -t147 * t129 - t150 * t179;
t1 = [(pkin(1) ^ 2 * t138 + t126 ^ 2 + t127 ^ 2) * MDP(7) + t110 ^ 2 * MDP(15) + MDP(1) + t124 ^ 2 * MDP(12) + (-0.2e1 * MDP(10) * t124 + MDP(8) * t115) * t115 + (-0.2e1 * MDP(23) * t106 + t176) * t107 + (0.2e1 * MDP(11) * t124 - 0.2e1 * MDP(9) * t115 + t169 + 0.2e1 * t170) * t114 + (-0.2e1 * MDP(18) * t114 - 0.2e1 * t171 + t172 + 0.2e1 * t173 - 0.2e1 * t174) * t109 + 0.2e1 * (-t127 * t183 - t142 * t191) * MDP(5) + 0.2e1 * (t126 * t183 + t145 * t191) * MDP(4) + (t106 * t96 + t109 * t94) * t193 + (t107 * t96 - t109 * t95) * t192 + 0.2e1 * (t102 * t109 + t114 * t98) * MDP(20) + 0.2e1 * (t102 * t110 - t114 * t99) * MDP(21) + 0.2e1 * (t105 * t124 + t108 * t115) * MDP(14) + 0.2e1 * (-t104 * t124 + t108 * t114) * MDP(13) + 0.2e1 * (-t126 * t142 + t127 * t145) * MDP(6) * t144; (t146 * t114 - t124 * t179) * MDP(13) + (t115 * t146 + t124 * t180) * MDP(14) + (-t109 * t179 - t128 * t114) * MDP(20) + (-t110 * t179 - t129 * t114) * MDP(21) + (t106 * t128 + t109 * t117) * MDP(27) + (t107 * t128 - t109 * t118) * MDP(28) + (-MDP(4) * t145 + MDP(5) * t142 - MDP(7) * pkin(1)) * t144; MDP(7); -t110 * t186 + t115 * MDP(10) - t114 * MDP(11) - t124 * MDP(12) + t104 * MDP(13) - t105 * MDP(14) + (t158 - t187) * t109 + (-t102 * MDP(20) + t162 * t114 - t153 + t171) * t151 + (t110 * MDP(15) + t102 * MDP(21) + t107 * t175 + (-t106 * t150 - t107 * t147) * MDP(23) + (pkin(9) * t106 + t185) * MDP(27) + (pkin(9) * t107 + t184) * MDP(28) + t163 * t114 + t155 * t109) * t148; (-t117 * t151 + t147 * t182) * MDP(27) + (t118 * t151 + t150 * t182) * MDP(28) + (-t149 * MDP(14) + (MDP(20) * t151 - MDP(21) * t148 + MDP(13)) * t152) * t143; MDP(12) + (t168 + 0.2e1 * t187) * t151 + (MDP(22) * t141 + MDP(15) - 0.2e1 * t165) * t140 + (-t121 * t151 + t140 * t190) * t193 + (t122 * t151 + t140 * t189) * t192 + 0.2e1 * (-t151 * t155 - t186) * t148; t170 + t169 + t98 * MDP(20) - t99 * MDP(21) + t147 * t176 + (-t106 * t147 + t107 * t150) * MDP(23) + (-pkin(4) * t106 - t184) * MDP(27) + (-pkin(4) * t107 + t185) * MDP(28) + (-MDP(18) + t154) * t109; -t129 * MDP(21) + (-MDP(20) - t157) * t128; (-t154 + t162) * t151 + (t147 * t175 + (-t139 + t141) * MDP(23) + (-pkin(4) * t147 - t189) * MDP(27) + (-pkin(4) * t150 + t190) * MDP(28) + t163) * t148; t139 * MDP(22) + 0.2e1 * pkin(4) * t157 + MDP(19) + 0.2e1 * t165; t153; MDP(27) * t117 - MDP(28) * t118; t160 * t148 + t158 - t168; t154; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
