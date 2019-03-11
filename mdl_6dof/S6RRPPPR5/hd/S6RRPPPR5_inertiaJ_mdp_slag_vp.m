% Calculate joint inertia matrix for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:46
% EndTime: 2019-03-09 08:23:48
% DurationCPUTime: 0.64s
% Computational Cost: add. (527->176), mult. (970->225), div. (0->0), fcn. (856->6), ass. (0->70)
t171 = pkin(4) + qJ(3);
t128 = sin(pkin(9));
t129 = cos(pkin(9));
t170 = t128 ^ 2 + t129 ^ 2;
t164 = pkin(3) + qJ(5);
t169 = t164 * t129;
t133 = cos(qJ(2));
t131 = sin(qJ(2));
t160 = t129 * t131;
t168 = pkin(4) * t160 + t133 * qJ(5);
t167 = -t164 * t128 - pkin(7);
t166 = 0.2e1 * t131;
t165 = pkin(7) * t133;
t163 = -pkin(5) - qJ(4);
t162 = pkin(2) * MDP(14);
t161 = t128 * t131;
t132 = cos(qJ(6));
t159 = t129 * t132;
t130 = sin(qJ(6));
t99 = t130 * t161 - t131 * t159;
t158 = t99 * MDP(26);
t109 = -pkin(2) * t133 - qJ(3) * t131 - pkin(1);
t96 = t128 * t109 + t129 * t165;
t157 = t170 * qJ(3) ^ 2;
t110 = t171 * t128;
t111 = t171 * t129;
t156 = MDP(27) * t133;
t105 = t128 * t132 + t129 * t130;
t100 = t105 * t131;
t155 = t100 * MDP(25);
t146 = qJ(4) * t128 + pkin(2);
t101 = t146 + t169;
t154 = t101 * MDP(22);
t153 = t105 * MDP(28);
t106 = t128 * t130 - t159;
t152 = t106 * MDP(23);
t107 = -pkin(3) * t129 - t146;
t151 = t107 * MDP(18);
t150 = MDP(12) - MDP(17);
t149 = MDP(15) - MDP(20);
t148 = MDP(18) + MDP(22);
t147 = MDP(21) - MDP(16);
t145 = MDP(11) + t147;
t144 = -MDP(19) + t150;
t114 = t128 * t165;
t95 = t109 * t129 - t114;
t93 = qJ(4) * t133 - t96;
t124 = t133 * pkin(3);
t94 = t124 - t95;
t143 = t128 * t94 - t129 * t93;
t142 = -t128 * t95 + t129 * t96;
t84 = t114 + t124 + (pkin(8) * t131 - t109) * t129 + t168;
t85 = t163 * t133 + (-pkin(4) - pkin(8)) * t161 + t96;
t141 = (-t130 * t84 + t132 * t85) * MDP(28) - (t130 * t85 + t132 * t84) * MDP(29);
t140 = t99 * MDP(28) + t100 * MDP(29);
t139 = MDP(28) * t132 - t130 * MDP(29);
t138 = -MDP(28) * t130 - MDP(29) * t132;
t137 = pkin(2) * MDP(11) + t107 * MDP(16) + t101 * MDP(21);
t136 = -pkin(2) * MDP(12) - t107 * MDP(17) + t101 * MDP(19);
t102 = pkin(8) * t128 + t110;
t103 = pkin(8) * t129 + t111;
t135 = t106 * MDP(25) + t105 * MDP(26) + (-t102 * t130 + t103 * t132) * MDP(28) - (t102 * t132 + t103 * t130) * MDP(29);
t113 = qJ(4) * t160;
t98 = -t113 + (pkin(3) * t128 + pkin(7)) * t131;
t97 = t163 * t128 - pkin(2) - t169;
t91 = -t167 * t131 - t113;
t90 = -pkin(4) * t161 - t93;
t89 = -t113 + (-pkin(5) * t129 - t167) * t131;
t86 = t94 + t168;
t1 = [(t95 ^ 2 + t96 ^ 2) * MDP(14) + (t93 ^ 2 + t94 ^ 2 + t98 ^ 2) * MDP(18) + (t86 ^ 2 + t90 ^ 2 + t91 ^ 2) * MDP(22) + MDP(1) + (MDP(23) * t100 - 0.2e1 * MDP(24) * t99) * t100 + (MDP(5) * t166 + 0.2e1 * pkin(1) * MDP(9) - 0.2e1 * t155 + t156 + 0.2e1 * t158) * t133 + 0.2e1 * t140 * t89 + 0.2e1 * (-t95 * MDP(11) + t96 * MDP(12) - t94 * MDP(16) + t93 * MDP(17) - t90 * MDP(19) + t86 * MDP(21) - t141) * t133 + ((-t95 * MDP(13) + t94 * MDP(15) - t98 * MDP(17) - t91 * MDP(19) - t86 * MDP(20)) * t129 + (-t96 * MDP(13) + t93 * MDP(15) - t98 * MDP(16) + t90 * MDP(20) + t91 * MDP(21)) * t128) * t166 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(4) + pkin(7) ^ 2 * MDP(14) + 0.2e1 * (MDP(11) * t128 + MDP(12) * t129) * pkin(7)) * t131) * t131; t142 * MDP(13) + t143 * MDP(15) + (-t128 * t86 - t129 * t90) * MDP(20) + (t110 * t86 + t111 * t90) * MDP(22) + t100 * t152 + (t100 * t105 - t106 * t99) * MDP(24) + (-t105 * t89 + t97 * t99) * MDP(28) + (t100 * t97 + t106 * t89) * MDP(29) + (-t128 * MDP(19) - t129 * MDP(21) - t154) * t91 + (t129 * MDP(16) - t128 * MDP(17) + t151) * t98 + (MDP(14) * t142 + MDP(18) * t143) * qJ(3) + (-pkin(7) * MDP(10) - t111 * MDP(19) + t110 * MDP(21) + MDP(7) + (t150 * t129 + (MDP(11) - MDP(16)) * t128) * qJ(3) - t135) * t133 + (MDP(6) + (-t110 * MDP(20) + t136) * t129 + (t111 * MDP(20) - t137) * t128 + (-t129 * MDP(11) + t128 * MDP(12) - MDP(9) - t162) * pkin(7)) * t131; MDP(8) + (pkin(2) ^ 2 + t157) * MDP(14) + (t107 ^ 2 + t157) * MDP(18) + (t101 ^ 2 + t110 ^ 2 + t111 ^ 2) * MDP(22) - 0.2e1 * t97 * t153 + (0.2e1 * t105 * MDP(24) + 0.2e1 * t97 * MDP(29) + t152) * t106 + 0.2e1 * (-t110 * t128 - t111 * t129) * MDP(20) + 0.2e1 * (MDP(13) + MDP(15)) * t170 * qJ(3) + 0.2e1 * t136 * t128 + 0.2e1 * t137 * t129; t98 * MDP(18) + MDP(22) * t91 + (pkin(7) * MDP(14) + t128 * t145 + t144 * t129) * t131 + t140; t106 * MDP(29) + t128 * t144 - t145 * t129 + t151 - t153 - t154 - t162; MDP(14) + t148; t94 * MDP(18) + t86 * MDP(22) + t149 * t160 + (-t138 + t147) * t133; t110 * MDP(22) + (MDP(18) * qJ(3) + t149) * t128; 0; t148; MDP(20) * t161 + t90 * MDP(22) + (-MDP(19) - t139) * t133; -t129 * MDP(20) + t111 * MDP(22); 0; 0; MDP(22); t141 + t155 - t156 - t158; t135; 0; t138; t139; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
