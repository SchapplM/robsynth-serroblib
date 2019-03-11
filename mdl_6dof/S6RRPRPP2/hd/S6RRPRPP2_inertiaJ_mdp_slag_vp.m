% Calculate joint inertia matrix for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:31
% EndTime: 2019-03-09 09:52:32
% DurationCPUTime: 0.58s
% Computational Cost: add. (964->184), mult. (1679->239), div. (0->0), fcn. (1744->6), ass. (0->66)
t122 = sin(pkin(9));
t123 = cos(pkin(9));
t125 = sin(qJ(2));
t160 = cos(qJ(2));
t107 = t122 * t160 + t123 * t125;
t164 = 0.2e1 * t107;
t163 = qJ(3) + pkin(7);
t106 = t122 * t125 - t123 * t160;
t124 = sin(qJ(4));
t126 = cos(qJ(4));
t115 = -t160 * pkin(2) - pkin(1);
t90 = t106 * pkin(3) - t107 * pkin(8) + t115;
t109 = t163 * t160;
t139 = t163 * t125;
t96 = t123 * t109 - t122 * t139;
t85 = -t124 * t96 + t126 * t90;
t83 = -t106 * pkin(4) - t85;
t86 = t124 * t90 + t126 * t96;
t162 = t85 * MDP(18) - t86 * MDP(19) - MDP(20) * t83;
t127 = pkin(4) + pkin(5);
t155 = qJ(5) * t126;
t161 = t124 * t127 - t155;
t159 = pkin(4) * MDP(20);
t158 = pkin(4) * MDP(23);
t114 = -pkin(2) * t123 - pkin(3);
t118 = t124 * qJ(5);
t150 = t126 * pkin(4) + t118;
t101 = t114 - t150;
t97 = pkin(5) * t126 - t101;
t156 = MDP(25) * t97;
t113 = pkin(2) * t122 + pkin(8);
t154 = t106 * t113;
t153 = t107 * t126;
t152 = -qJ(6) + t113;
t151 = t126 * MDP(24) + t124 * MDP(25);
t120 = t124 ^ 2;
t121 = t126 ^ 2;
t149 = t120 + t121;
t148 = MDP(14) * t126;
t147 = 0.2e1 * t160;
t146 = MDP(18) + MDP(20);
t145 = -MDP(19) + MDP(22);
t144 = MDP(20) + MDP(24);
t143 = MDP(21) - MDP(26);
t142 = MDP(22) + MDP(25);
t141 = MDP(23) + MDP(27);
t100 = t106 * qJ(5);
t140 = 0.2e1 * t100 + t86;
t82 = t100 + t86;
t94 = t109 * t122 + t123 * t139;
t138 = MDP(27) * t127 + MDP(20);
t137 = MDP(23) * t113 + MDP(21);
t136 = pkin(4) * t124 - t155;
t103 = t152 * t124;
t104 = t152 * t126;
t134 = -t103 * MDP(24) + t104 * MDP(25);
t84 = -t161 * t107 - t94;
t87 = t136 * t107 + t94;
t133 = -t94 * MDP(18) - t87 * MDP(20) + t84 * MDP(24);
t80 = -pkin(5) * t106 - qJ(6) * t153 + t83;
t132 = t94 * MDP(19) - t87 * MDP(22) + t84 * MDP(25) - t80 * MDP(26);
t131 = -MDP(18) * t114 - MDP(20) * t101 + MDP(24) * t97 - MDP(26) * t104;
t129 = qJ(5) ^ 2;
t98 = t124 * t107 * qJ(6);
t81 = t98 + t82;
t1 = [MDP(1) + pkin(1) * MDP(9) * t147 + (t115 ^ 2 + t94 ^ 2 + t96 ^ 2) * MDP(12) + (t82 ^ 2 + t83 ^ 2 + t87 ^ 2) * MDP(23) + (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) * MDP(27) + (MDP(13) * t121 - 0.2e1 * t124 * t148) * t107 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t125 + MDP(5) * t147) * t125 + (t106 * MDP(17) + (t126 * MDP(15) - t124 * MDP(16)) * t164) * t106 + 0.2e1 * (-MDP(11) * t96 + MDP(22) * t82 - MDP(24) * t80 + MDP(25) * t81 + t162) * t106 + (t94 * MDP(11) + (t83 * MDP(21) + t132) * t126 + (-t82 * MDP(21) + t81 * MDP(26) - t133) * t124) * t164; t125 * MDP(6) + t160 * MDP(7) + t87 * t101 * MDP(23) + (t103 * t80 + t104 * t81 + t84 * t97) * MDP(27) + (-t120 + t121) * MDP(14) * t107 + t134 * t106 + (-t160 * MDP(10) - t125 * MDP(9)) * pkin(7) + ((-t106 * t122 - t107 * t123) * MDP(11) + (t122 * t96 - t123 * t94) * MDP(12)) * pkin(2) + (t106 * MDP(16) + (t107 * t114 - t154) * MDP(19) + (-t101 * t107 + t154) * MDP(22) + t107 * t156 + (-t103 * t107 - t81) * MDP(26) + t137 * t82 + t133) * t126 + (t137 * t83 + (-t146 * t113 + MDP(15)) * t106 + (MDP(13) * t126 - t131) * t107 + t132) * t124; MDP(8) + t120 * MDP(13) + (t149 * t113 ^ 2 + t101 ^ 2) * MDP(23) + (t103 ^ 2 + t104 ^ 2 + t97 ^ 2) * MDP(27) + (t122 ^ 2 + t123 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t149 * MDP(21) * t113 + 0.2e1 * t131 * t126 + 0.2e1 * (MDP(19) * t114 - MDP(22) * t101 - MDP(26) * t103 + t148 + t156) * t124; t115 * MDP(12) + (t124 * t82 - t126 * t83) * MDP(23) + (t124 * t81 - t126 * t80) * MDP(27) - t143 * t149 * t107 + ((MDP(18) + t144) * t126 + (-MDP(19) + t142) * t124) * t106; (-t103 * t126 + t104 * t124) * MDP(27); t141 * t149 + MDP(12); t140 * MDP(22) + (-pkin(4) * t83 + qJ(5) * t82) * MDP(23) - t83 * MDP(24) + (t98 + t140) * MDP(25) + (qJ(5) * t81 - t127 * t80) * MDP(27) + (MDP(17) + t159 + (pkin(5) + t127) * MDP(24)) * t106 + ((-t143 * qJ(5) - MDP(16)) * t124 + (-MDP(21) * pkin(4) + MDP(24) * qJ(6) + MDP(26) * t127 + MDP(15)) * t126) * t107 + t162; t124 * MDP(15) + t126 * MDP(16) - t136 * MDP(21) + t161 * MDP(26) + (qJ(5) * t104 - t103 * t127) * MDP(27) + ((MDP(23) * qJ(5) + t145) * t126 + (-t146 - t158) * t124) * t113 + t134; t150 * MDP(23) + t118 * MDP(27) + (MDP(18) + t138) * t126 + t145 * t124 + t151; MDP(17) + 0.2e1 * t159 + (pkin(4) ^ 2 + t129) * MDP(23) + 0.2e1 * t127 * MDP(24) + (t127 ^ 2 + t129) * MDP(27) + 0.2e1 * t142 * qJ(5); MDP(23) * t83 + MDP(27) * t80 - t144 * t106 + t143 * t153; MDP(27) * t103 + (-MDP(26) + t137) * t124; -t141 * t126; -MDP(24) - t138 - t158; t141; MDP(27) * t84 + (-t124 * MDP(24) + t126 * MDP(25)) * t107; MDP(27) * t97 + t151; 0; 0; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
