% Calculate joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR9_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:59
% EndTime: 2021-01-15 23:24:03
% DurationCPUTime: 0.63s
% Computational Cost: add. (738->155), mult. (1485->234), div. (0->0), fcn. (1525->8), ass. (0->71)
t138 = sin(qJ(3));
t166 = qJ(4) + pkin(7);
t123 = t166 * t138;
t141 = cos(qJ(3));
t124 = t166 * t141;
t135 = sin(pkin(9));
t136 = cos(pkin(9));
t106 = -t136 * t123 - t124 * t135;
t107 = -t123 * t135 + t124 * t136;
t117 = t135 * t138 - t136 * t141;
t118 = t135 * t141 + t136 * t138;
t137 = sin(qJ(5));
t140 = cos(qJ(5));
t100 = t117 * t140 + t118 * t137;
t101 = -t117 * t137 + t118 * t140;
t94 = -pkin(8) * t118 + t106;
t95 = -pkin(8) * t117 + t107;
t149 = t101 * MDP(24) - t100 * MDP(25) + (-t137 * t95 + t140 * t94) * MDP(27) - (t137 * t94 + t140 * t95) * MDP(28);
t176 = t149 + t138 * MDP(13) + t141 * MDP(14) + t106 * MDP(18) - t107 * MDP(19) - (MDP(16) * t138 + MDP(17) * t141) * pkin(7);
t139 = sin(qJ(2));
t112 = t118 * t139;
t162 = t139 * t141;
t164 = t138 * t139;
t113 = -t135 * t164 + t136 * t162;
t92 = t112 * t140 + t113 * t137;
t93 = -t112 * t137 + t113 * t140;
t165 = t93 * MDP(24) - t92 * MDP(25);
t142 = cos(qJ(2));
t122 = -pkin(2) * t142 - pkin(7) * t139 - pkin(1);
t119 = t141 * t122;
t169 = pkin(6) * t138;
t102 = -qJ(4) * t162 + t119 + (-pkin(3) - t169) * t142;
t167 = pkin(6) * t142;
t151 = t141 * t167;
t105 = t151 + (-qJ(4) * t139 + t122) * t138;
t88 = t102 * t136 - t105 * t135;
t84 = -pkin(4) * t142 - pkin(8) * t113 + t88;
t89 = t102 * t135 + t136 * t105;
t87 = -pkin(8) * t112 + t89;
t79 = -t137 * t87 + t140 * t84;
t80 = t137 * t84 + t140 * t87;
t175 = t79 * MDP(27) - t80 * MDP(28) + t165;
t173 = 2 * MDP(20);
t172 = -2 * MDP(23);
t171 = 0.2e1 * MDP(28);
t170 = pkin(3) * t135;
t168 = pkin(6) * t141;
t163 = t138 * t141;
t121 = pkin(3) * t164 + t139 * pkin(6);
t159 = MDP(22) * t101;
t158 = t100 * MDP(27);
t128 = pkin(3) * t136 + pkin(4);
t157 = (t128 * t140 - t137 * t170) * MDP(27);
t156 = (t128 * t137 + t140 * t170) * MDP(28);
t155 = t117 * MDP(18);
t154 = t118 * MDP(19);
t129 = -pkin(3) * t141 - pkin(2);
t153 = t129 * MDP(21);
t152 = MDP(15) + MDP(26);
t150 = MDP(12) * t163;
t148 = MDP(13) * t141 - MDP(14) * t138;
t146 = MDP(18) * t136 - MDP(19) * t135;
t145 = MDP(26) - t156 + t157;
t133 = t141 ^ 2;
t132 = t139 ^ 2;
t131 = t138 ^ 2;
t111 = pkin(4) * t117 + t129;
t110 = t122 * t138 + t151;
t109 = -t138 * t167 + t119;
t103 = pkin(4) * t112 + t121;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t139 * MDP(10) + (t121 ^ 2 + t88 ^ 2 + t89 ^ 2) * MDP(21) + (MDP(22) * t93 + t172 * t92) * t93 + t152 * t142 ^ 2 + (MDP(11) * t133 + MDP(4) - 0.2e1 * t150) * t132 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t148) * t139 - t165) * t142 + 0.2e1 * (-t109 * t142 + t132 * t169) * MDP(16) + 0.2e1 * (t110 * t142 + t132 * t168) * MDP(17) + 0.2e1 * (t112 * t121 - t142 * t88) * MDP(18) + 0.2e1 * (t113 * t121 + t142 * t89) * MDP(19) + (-t112 * t89 - t113 * t88) * t173 + 0.2e1 * (t103 * t92 - t142 * t79) * MDP(27) + (t103 * t93 + t142 * t80) * t171; (t112 * t129 + t117 * t121) * MDP(18) + (t113 * t129 + t118 * t121) * MDP(19) + (-t106 * t113 - t107 * t112 - t117 * t89 - t118 * t88) * MDP(20) + (t106 * t88 + t107 * t89 + t121 * t129) * MDP(21) + t93 * t159 + (-t100 * t93 - t101 * t92) * MDP(23) + (t100 * t103 + t111 * t92) * MDP(27) + (t101 * t103 + t111 * t93) * MDP(28) + (-pkin(6) * MDP(10) + MDP(7) - t176) * t142 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t163 + (-t131 + t133) * MDP(12) + (-pkin(2) * t138 - t168) * MDP(16) + (-pkin(2) * t141 + t169) * MDP(17)) * t139; MDP(8) + t131 * MDP(11) + 0.2e1 * t150 + (-t106 * t118 - t107 * t117) * t173 + (t106 ^ 2 + t107 ^ 2) * MDP(21) + 0.2e1 * t111 * t158 + (t153 + 0.2e1 * t154 + 0.2e1 * t155) * t129 + 0.2e1 * (MDP(16) * t141 - MDP(17) * t138) * pkin(2) + (t100 * t172 + t111 * t171 + t159) * t101; t109 * MDP(16) - t110 * MDP(17) + t88 * MDP(18) - t89 * MDP(19) + (-MDP(15) - t145) * t142 + t148 * t139 + ((-t112 * t135 - t113 * t136) * MDP(20) + (t135 * t89 + t136 * t88) * MDP(21) - t146 * t142) * pkin(3) + t175; ((-t117 * t135 - t118 * t136) * MDP(20) + (t106 * t136 + t107 * t135) * MDP(21)) * pkin(3) + t176; 0.2e1 * t157 - 0.2e1 * t156 + t152 + (0.2e1 * t146 + (t135 ^ 2 + t136 ^ 2) * MDP(21) * pkin(3)) * pkin(3); MDP(18) * t112 + MDP(19) * t113 + MDP(21) * t121 + MDP(27) * t92 + MDP(28) * t93; MDP(28) * t101 + t153 + t154 + t155 + t158; 0; MDP(21); -t142 * MDP(26) + t175; t149; t145; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
