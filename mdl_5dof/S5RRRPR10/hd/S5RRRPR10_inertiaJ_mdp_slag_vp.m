% Calculate joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR10_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:20
% EndTime: 2021-01-15 23:41:24
% DurationCPUTime: 0.77s
% Computational Cost: add. (1104->186), mult. (2498->288), div. (0->0), fcn. (2712->10), ass. (0->85)
t140 = sin(pkin(5));
t188 = 0.2e1 * t140;
t147 = cos(qJ(3));
t179 = qJ(4) + pkin(8);
t129 = t179 * t147;
t139 = sin(pkin(10));
t141 = cos(pkin(10));
t144 = sin(qJ(3));
t161 = t179 * t144;
t115 = t141 * t129 - t139 * t161;
t187 = t144 * MDP(13) + t147 * MDP(14) - t115 * MDP(19) - (t144 * MDP(16) + t147 * MDP(17)) * pkin(8);
t186 = 0.2e1 * MDP(16);
t185 = 2 * MDP(18);
t184 = 2 * MDP(20);
t183 = 2 * MDP(27);
t182 = 2 * MDP(28);
t145 = sin(qJ(2));
t181 = pkin(1) * t145;
t148 = cos(qJ(2));
t180 = pkin(1) * t148;
t142 = cos(pkin(5));
t176 = t140 * t148;
t163 = pkin(7) * t176;
t118 = t163 + (pkin(8) + t181) * t142;
t119 = (-pkin(2) * t148 - pkin(8) * t145 - pkin(1)) * t140;
t107 = -t144 * t118 + t147 * t119;
t177 = t140 * t145;
t122 = t142 * t144 + t147 * t177;
t100 = -pkin(3) * t176 - t122 * qJ(4) + t107;
t108 = t147 * t118 + t144 * t119;
t121 = -t142 * t147 + t144 * t177;
t104 = -t121 * qJ(4) + t108;
t96 = t139 * t100 + t141 * t104;
t113 = t139 * t129 + t141 * t161;
t127 = t139 * t147 + t141 * t144;
t178 = t113 * t127;
t175 = t142 * MDP(8);
t174 = t145 * MDP(6);
t110 = -t139 * t121 + t141 * t122;
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t105 = t143 * t110 + t146 * t176;
t173 = t105 * MDP(23);
t172 = t105 * MDP(25);
t106 = t146 * t110 - t143 * t176;
t171 = t106 * MDP(22);
t109 = t141 * t121 + t139 * t122;
t170 = t109 * MDP(24);
t169 = t109 * MDP(26);
t168 = t122 * MDP(11);
t126 = t139 * t144 - t141 * t147;
t167 = t126 * MDP(26);
t166 = t127 * MDP(19);
t165 = t143 * MDP(22);
t164 = t148 * MDP(15);
t135 = -t147 * pkin(3) - pkin(2);
t162 = t143 * t146 * MDP(23);
t95 = t141 * t100 - t139 * t104;
t160 = t122 * MDP(13) - t121 * MDP(14);
t158 = t141 * MDP(18) - t139 * MDP(19);
t130 = pkin(7) * t177;
t117 = t130 + (-pkin(2) - t180) * t142;
t111 = t121 * pkin(3) + t117;
t157 = t110 * MDP(19) + t111 * MDP(21);
t156 = t146 * MDP(27) - t143 * MDP(28);
t155 = MDP(27) * t143 + MDP(28) * t146;
t154 = MDP(18) + t156;
t153 = (MDP(24) * t146 - MDP(25) * t143) * t127;
t151 = t143 * MDP(24) + t146 * MDP(25) - t155 * (t139 * pkin(3) + pkin(9));
t94 = -pkin(9) * t176 + t96;
t97 = t109 * pkin(4) - t110 * pkin(9) + t111;
t91 = -t143 * t94 + t146 * t97;
t92 = t143 * t97 + t146 * t94;
t150 = t106 * MDP(24) + t91 * MDP(27) - t92 * MDP(28) + t169 - t172;
t138 = t146 ^ 2;
t137 = t143 ^ 2;
t136 = t140 ^ 2;
t134 = -t141 * pkin(3) - pkin(4);
t124 = t142 * t181 + t163;
t123 = t142 * t180 - t130;
t112 = t126 * pkin(4) - t127 * pkin(9) + t135;
t102 = t143 * t112 + t146 * t115;
t101 = t146 * t112 - t143 * t115;
t93 = pkin(4) * t176 - t95;
t1 = [(t111 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(21) + MDP(1) + t136 * t145 ^ 2 * MDP(4) + (t174 * t188 + t175) * t142 + (-0.2e1 * t121 * MDP(12) + t168) * t122 + (t169 - 0.2e1 * t172) * t109 + (0.2e1 * t170 + t171 - 0.2e1 * t173) * t106 + ((0.2e1 * MDP(5) * t145 + t164) * t136 + (MDP(7) * t142 - t160) * t188) * t148 + (t93 * t105 + t91 * t109) * t183 + (t93 * t106 - t92 * t109) * t182 + (-t96 * t109 - t95 * t110) * t184 + 0.2e1 * (-t124 * t142 - t136 * t181) * MDP(10) + (t111 * t109 - t95 * t176) * t185 + 0.2e1 * (t111 * t110 + t96 * t176) * MDP(19) + (-t107 * t176 + t117 * t121) * t186 + 0.2e1 * (t108 * t176 + t117 * t122) * MDP(17) + 0.2e1 * (t123 * t142 + t136 * t180) * MDP(9); t175 + t123 * MDP(9) - t124 * MDP(10) + t144 * t168 + (-t144 * t121 + t122 * t147) * MDP(12) + (-pkin(2) * t121 - t117 * t147) * MDP(16) + (-pkin(2) * t122 + t117 * t144) * MDP(17) + (-t115 * t109 + t113 * t110) * MDP(20) + (-t95 * t113 + t96 * t115) * MDP(21) + (t101 * t109 + t113 * t105) * MDP(27) + (-t102 * t109 + t113 * t106) * MDP(28) + (t109 * MDP(18) + t157) * t135 + (t111 * MDP(18) - t96 * MDP(20) + t150) * t126 + (t111 * MDP(19) - t95 * MDP(20) + (-t106 * MDP(23) - t109 * MDP(25) + t93 * MDP(27)) * t143 + (t93 * MDP(28) + t170 + t171 - t173) * t146) * t127 + (t174 + (t113 * MDP(18) + MDP(7) - t187) * t148) * t140; MDP(8) + pkin(2) * t147 * t186 + 0.2e1 * t135 * t166 + (t113 ^ 2 + t115 ^ 2 + t135 ^ 2) * MDP(21) + (t138 * MDP(22) - 0.2e1 * t162) * t127 ^ 2 + (MDP(11) * t144 + 0.2e1 * t147 * MDP(12) - 0.2e1 * pkin(2) * MDP(17)) * t144 + (t135 * t185 + 0.2e1 * t153 + t167) * t126 + (-t115 * t126 + t178) * t184 + (t101 * t126 + t143 * t178) * t183 + (-t102 * t126 + t146 * t178) * t182; -t140 * t164 + t107 * MDP(16) - t108 * MDP(17) + t95 * MDP(18) - t96 * MDP(19) + t106 * t165 + (-t143 * t105 + t106 * t146) * MDP(23) + (t134 * t105 - t93 * t146) * MDP(27) + (t134 * t106 + t93 * t143) * MDP(28) + t151 * t109 + ((-t109 * t139 - t110 * t141) * MDP(20) + (t139 * t96 + t141 * t95) * MDP(21) - t158 * t176) * pkin(3) + t160; -t154 * t113 + t151 * t126 + (t146 * t165 + (-t137 + t138) * MDP(23) + t155 * t134) * t127 + ((-t126 * t139 - t127 * t141) * MDP(20) + (-t113 * t141 + t115 * t139) * MDP(21)) * pkin(3) + t187; t137 * MDP(22) - 0.2e1 * t156 * t134 + MDP(15) + 0.2e1 * t162 + (0.2e1 * t158 + (t139 ^ 2 + t141 ^ 2) * MDP(21) * pkin(3)) * pkin(3); t154 * t109 + t157; t135 * MDP(21) + t154 * t126 + t166; 0; MDP(21); t150; t101 * MDP(27) - t102 * MDP(28) + t153 + t167; t151; t156; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
