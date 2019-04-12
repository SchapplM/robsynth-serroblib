% Calculate joint inertia matrix for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR14V3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14V3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_inertiaJ_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR14V3_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:09
% EndTime: 2019-04-12 15:10:12
% DurationCPUTime: 0.85s
% Computational Cost: add. (333->163), mult. (1015->239), div. (0->0), fcn. (1011->8), ass. (0->76)
t112 = cos(qJ(4));
t106 = sin(qJ(6));
t110 = cos(qJ(6));
t111 = cos(qJ(5));
t157 = t111 * t112;
t108 = sin(qJ(4));
t160 = t108 * t110;
t163 = t106 * t108;
t159 = t108 * t111;
t93 = t106 * t159 + t110 * t112;
t96 = -t106 * t112 + t110 * t159;
t117 = ((-t106 * t157 + t160) * MDP(34) - (t110 * t157 + t163) * MDP(35)) * qJ(3) + t96 * MDP(31) - t93 * MDP(32);
t115 = t112 * MDP(25) + t117;
t174 = 0.2e1 * qJ(3);
t173 = 0.2e1 * t108;
t172 = t93 * MDP(34);
t107 = sin(qJ(5));
t100 = t107 ^ 2;
t124 = MDP(34) * t106 + MDP(35) * t110;
t171 = t124 * t100 - MDP(21);
t168 = -2 * MDP(30);
t113 = cos(qJ(2));
t109 = sin(qJ(2));
t158 = t109 * t112;
t94 = t107 * t158 + t111 * t113;
t167 = t107 * t94;
t166 = qJ(3) * t109;
t165 = qJ(3) * t112;
t101 = t108 ^ 2;
t164 = t101 * t111;
t162 = t107 * t111;
t161 = t108 * t109;
t137 = t109 * t160;
t136 = t109 * t157;
t95 = -t107 * t113 + t136;
t87 = t106 * t95 - t137;
t156 = t87 * MDP(32);
t88 = t106 * t161 + t110 * t95;
t155 = t88 * MDP(31);
t154 = t94 * MDP(33);
t153 = t95 * MDP(22);
t152 = t95 * MDP(23);
t151 = t95 * MDP(24);
t150 = t96 * MDP(29);
t105 = t112 ^ 2;
t149 = t101 + t105;
t148 = MDP(28) * t111;
t147 = t106 * MDP(29);
t146 = t107 * MDP(28);
t145 = t108 * MDP(24);
t144 = t110 * MDP(29);
t143 = t112 * MDP(17);
t141 = t113 * MDP(19);
t140 = qJ(3) ^ 2 * MDP(14);
t135 = MDP(30) * t106 * t110;
t134 = MDP(27) * t158;
t133 = t108 * MDP(21) - MDP(11);
t132 = t112 * MDP(21) + MDP(12);
t131 = -t111 * MDP(24) + MDP(16);
t129 = t87 * MDP(34) + t88 * MDP(35);
t90 = (-t106 * t158 + t111 * t137) * qJ(3);
t128 = t90 * MDP(35) + t166 * t172;
t126 = t110 * MDP(31) - t106 * MDP(32);
t125 = t106 * MDP(31) + t110 * MDP(32);
t123 = t111 * MDP(27) + MDP(20) - t146;
t122 = -MDP(23) + t126;
t121 = -MDP(25) + t125;
t120 = MDP(34) * t110 - MDP(35) * t106 + MDP(27);
t119 = t154 + t155 - t156;
t118 = t123 * t112;
t116 = t119 + t128;
t104 = t111 ^ 2;
t103 = t110 ^ 2;
t102 = t109 ^ 2;
t99 = t106 ^ 2;
t1 = [MDP(1) + (0.2e1 * t109 * t145 + t153) * t95 + (t88 * MDP(29) + t87 * t168) * t88 + (-0.2e1 * t108 * t112 * MDP(16) + MDP(13) * t174 + t105 * MDP(15) + t101 * MDP(26) + MDP(4) + t140) * t102 + ((t102 * t164 + t95 * t158) * MDP(28) + (t101 * t102 * MDP(27) - t129 * t161) * t107) * t174 + (t141 + (-0.2e1 * t143 + MDP(18) * t173 + (2 * MDP(5)) + (t112 * MDP(20) - t133) * t174) * t109) * t113 + (-0.2e1 * MDP(25) * t161 + t134 * t174 + 0.2e1 * t128 - 0.2e1 * t152 + t154 + 0.2e1 * t155 - 0.2e1 * t156) * t94; t88 * t150 + (-t87 * t96 - t88 * t93) * MDP(30) + (t129 * t107 * qJ(3) - t151) * t112 + t115 * t94 + (t105 * MDP(16) + MDP(6) + (-t107 * MDP(25) - t131) * t101) * t109 + (-t112 * MDP(18) + t132 * qJ(3) + MDP(7)) * t113 + (-t113 * MDP(17) + (-t94 * MDP(23) + t153) * t111 + (MDP(15) - MDP(26)) * t158 + (t113 * MDP(20) + t94 * MDP(27) + (t95 - t136) * MDP(28)) * qJ(3) + (-t152 - qJ(3) * t134 + (-t96 * t166 + t90) * MDP(35) + t119) * t107) * t108; t140 + t105 * MDP(26) + MDP(8) + (t93 * t168 + t150) * t96 + t131 * t112 * t173 + (t104 * MDP(22) + t100 * MDP(33) + MDP(15)) * t101 + (t149 * t148 + MDP(13)) * t174 + (-0.2e1 * MDP(23) * t164 + (t149 * MDP(27) + (t96 * MDP(35) + t172) * t112) * t174 + t115 * t173) * t107; (-t106 * t167 - t111 * t87) * MDP(34) + (-t110 * t167 - t111 * t88) * MDP(35) + (t123 * t108 + t132) * t109; (-t100 * t163 - t111 * t93) * MDP(34) + (-t100 * t160 - t111 * t96) * MDP(35) - t118 + t133; MDP(14); -t141 + (-t116 + t152) * t111 + (t153 + t88 * t144 + (-t106 * t88 - t110 * t87) * MDP(30) + t122 * t94) * t107 + (t143 + (t107 * MDP(24) + t111 * MDP(25) - MDP(18)) * t108 + (-t171 * t108 - t118) * qJ(3)) * t109; (t171 * qJ(3) + MDP(18)) * t112 - t115 * t111 + (-t112 * MDP(24) + t96 * t144 + (-t106 * t96 - t110 * t93) * MDP(30)) * t107 + (t104 * MDP(23) + MDP(17) + (MDP(22) - MDP(33)) * t162 + t122 * t100 - t123 * qJ(3)) * t108; 0; MDP(33) * t104 + MDP(19) - 0.2e1 * t122 * t162 + (MDP(29) * t103 + MDP(22) - 0.2e1 * t135) * t100; t151 + t88 * t147 + (-t106 * t87 + t110 * t88) * MDP(30) + t121 * t94 + (MDP(26) + (t120 * t107 + t148) * qJ(3)) * t161; -t112 * MDP(26) + t96 * t147 + (-t106 * t93 + t110 * t96) * MDP(30) + (-MDP(28) * t165 + t145) * t111 + (t121 * t108 - t120 * t165) * t107; t120 * t111 - t146; -t121 * t111 + (MDP(24) + t106 * t144 + (t103 - t99) * MDP(30)) * t107; MDP(29) * t99 + MDP(26) + 0.2e1 * t135; t116; t108 * t107 * MDP(33) + t117; -t124 * t107; -t111 * MDP(33) + t126 * t107; t125; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
