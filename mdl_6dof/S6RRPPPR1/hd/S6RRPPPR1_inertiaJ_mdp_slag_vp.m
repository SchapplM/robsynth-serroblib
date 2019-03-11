% Calculate joint inertia matrix for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:20
% EndTime: 2019-03-09 08:08:21
% DurationCPUTime: 0.60s
% Computational Cost: add. (838->146), mult. (1529->209), div. (0->0), fcn. (1645->8), ass. (0->72)
t151 = MDP(16) + MDP(20);
t133 = sin(pkin(10));
t135 = cos(pkin(10));
t160 = t133 ^ 2 + t135 ^ 2;
t182 = t151 * t160;
t181 = qJ(3) + pkin(7);
t180 = (MDP(15) + MDP(18)) * t160;
t137 = sin(qJ(6));
t139 = cos(qJ(6));
t179 = -t133 * t139 + t135 * t137;
t153 = MDP(14) - MDP(19);
t154 = MDP(13) + MDP(17);
t178 = t133 * t154 + t135 * t153;
t118 = t133 * t137 + t135 * t139;
t113 = t118 * MDP(26);
t162 = MDP(27) * t179 - t113;
t177 = -t133 * t153 + t135 * t154 - t162;
t171 = cos(qJ(2));
t124 = t181 * t171;
t134 = sin(pkin(9));
t136 = cos(pkin(9));
t138 = sin(qJ(2));
t150 = t181 * t138;
t104 = t124 * t134 + t136 * t150;
t176 = t104 ^ 2;
t129 = -pkin(2) * t136 - pkin(3);
t142 = qJ(5) * t133 - t129;
t172 = pkin(4) + pkin(5);
t108 = t135 * t172 + t142;
t175 = 0.2e1 * t108;
t174 = -0.2e1 * t135;
t173 = -2 * MDP(22);
t127 = pkin(2) * t134 + qJ(4);
t170 = -pkin(8) + t127;
t116 = t134 * t138 - t136 * t171;
t119 = t134 * t171 + t136 * t138;
t130 = -pkin(2) * t171 - pkin(1);
t100 = t116 * pkin(3) - t119 * qJ(4) + t130;
t106 = t136 * t124 - t134 * t150;
t93 = t133 * t100 + t135 * t106;
t169 = qJ(5) * t135;
t168 = t119 * t133;
t167 = t119 * t135;
t95 = t179 * t119;
t164 = t95 * MDP(24);
t96 = t118 * t119;
t163 = t96 * MDP(23);
t159 = MDP(25) * t116;
t110 = -pkin(4) * t135 - t142;
t158 = t110 * MDP(20);
t157 = t179 * MDP(21);
t156 = t129 * MDP(16);
t155 = 0.2e1 * t171;
t89 = t116 * qJ(5) + t93;
t102 = t133 * t106;
t92 = t100 * t135 - t102;
t90 = -pkin(4) * t116 - t92;
t149 = t133 * t90 + t135 * t89;
t148 = t133 * t89 - t135 * t90;
t147 = t133 * t93 + t135 * t92;
t146 = -t133 * t92 + t135 * t93;
t87 = t102 + (-pkin(8) * t119 - t100) * t135 - t172 * t116;
t88 = pkin(8) * t168 + t89;
t145 = (-t137 * t88 + t139 * t87) * MDP(26) - (t137 * t87 + t139 * t88) * MDP(27);
t144 = MDP(26) * t95 + MDP(27) * t96;
t143 = MDP(26) * t139 - t137 * MDP(27);
t111 = t170 * t133;
t112 = t170 * t135;
t141 = -t179 * MDP(23) - t118 * MDP(24) + (t111 * t139 - t112 * t137) * MDP(26) - (t111 * t137 + t112 * t139) * MDP(27);
t94 = (pkin(4) * t133 - t169) * t119 + t104;
t91 = (-t133 * t172 + t169) * t119 - t104;
t1 = [MDP(1) + pkin(1) * MDP(9) * t155 + (t106 ^ 2 + t130 ^ 2 + t176) * MDP(12) + (t92 ^ 2 + t93 ^ 2 + t176) * MDP(16) + (t89 ^ 2 + t90 ^ 2 + t94 ^ 2) * MDP(20) + (t96 * MDP(21) + t173 * t95) * t96 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t138 + MDP(5) * t155) * t138 + (t159 - 0.2e1 * t163 + 0.2e1 * t164) * t116 + 0.2e1 * t144 * t91 + 0.2e1 * (-t106 * MDP(11) + t92 * MDP(13) - t93 * MDP(14) - t90 * MDP(17) + t89 * MDP(19) - t145) * t116 + 0.2e1 * (-t147 * MDP(15) - t148 * MDP(18) + (MDP(17) * t133 - MDP(19) * t135) * t94 + (MDP(13) * t133 + MDP(14) * t135 + MDP(11)) * t104) * t119; t138 * MDP(6) + t171 * MDP(7) + (-t104 * t135 + t129 * t168) * MDP(13) + (t104 * t133 + t129 * t167) * MDP(14) + t146 * MDP(15) + t104 * t156 + (t110 * t168 - t135 * t94) * MDP(17) + t149 * MDP(18) + (-t110 * t167 - t133 * t94) * MDP(19) + t94 * t158 - t96 * t157 + (-t96 * t118 + t179 * t95) * MDP(22) + (t108 * t95 + t118 * t91) * MDP(26) + (t108 * t96 - t179 * t91) * MDP(27) + (t146 * MDP(16) + t149 * MDP(20)) * t127 + (-t178 * t127 - t141) * t116 + (-MDP(10) * t171 - t138 * MDP(9)) * pkin(7) + ((-t116 * t134 - t119 * t136) * MDP(11) + (-t104 * t136 + t106 * t134) * MDP(12)) * pkin(2); MDP(8) + t113 * t175 + (t134 ^ 2 + t136 ^ 2) * MDP(12) * pkin(2) ^ 2 + (MDP(13) * t174 + 0.2e1 * t133 * MDP(14) + t156) * t129 + (MDP(17) * t174 - 0.2e1 * t133 * MDP(19) + t158) * t110 - (MDP(27) * t175 + t118 * t173 - t157) * t179 + (t182 * t127 + 0.2e1 * t180) * t127; t130 * MDP(12) + t147 * MDP(16) + t148 * MDP(20) + t177 * t116 - t119 * t180; 0; MDP(12) + t182; t104 * MDP(16) + t94 * MDP(20) + t178 * t119 - t144; t156 + t158 - t177; 0; t151; MDP(18) * t167 + t90 * MDP(20) + (-MDP(17) - t143) * t116; (MDP(20) * t127 + MDP(18)) * t133; -t135 * MDP(20); 0; MDP(20); t145 - t159 + t163 - t164; t141; t162; 0; t143; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
