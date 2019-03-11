% Calculate joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:08
% EndTime: 2019-03-08 21:54:09
% DurationCPUTime: 0.47s
% Computational Cost: add. (683->124), mult. (1398->192), div. (0->0), fcn. (1658->12), ass. (0->73)
t131 = sin(qJ(6));
t135 = cos(qJ(6));
t143 = t135 * MDP(26) - t131 * MDP(27);
t170 = MDP(19) + t143;
t158 = t131 * MDP(23) + t135 * MDP(24);
t137 = cos(qJ(3));
t121 = -t137 * pkin(3) - pkin(2);
t127 = sin(pkin(12));
t129 = cos(pkin(12));
t133 = sin(qJ(3));
t146 = t127 * t133 - t129 * t137;
t105 = pkin(4) * t146 + t121;
t169 = 0.2e1 * t105;
t168 = t127 * pkin(3);
t166 = -qJ(4) - pkin(8);
t165 = MDP(13) * pkin(3);
t132 = sin(qJ(5));
t136 = cos(qJ(5));
t117 = t166 * t133;
t118 = t166 * t137;
t103 = t129 * t117 + t127 * t118;
t113 = t127 * t137 + t129 * t133;
t91 = -t113 * pkin(9) + t103;
t104 = t127 * t117 - t129 * t118;
t92 = -pkin(9) * t146 + t104;
t84 = t132 * t92 - t136 * t91;
t164 = t84 * t135;
t102 = t136 * t113 - t132 * t146;
t163 = t102 * t131;
t162 = t102 * t135;
t128 = sin(pkin(6));
t134 = sin(qJ(2));
t161 = t128 * t134;
t138 = cos(qJ(2));
t160 = t128 * t138;
t159 = t131 * t135;
t157 = MDP(10) * t137;
t101 = t132 * t113 + t136 * t146;
t156 = t101 * MDP(25);
t155 = t102 * MDP(20);
t120 = t129 * pkin(3) + pkin(4);
t108 = t136 * t120 - t132 * t168;
t154 = t108 * MDP(19);
t109 = -t132 * t120 - t136 * t168;
t153 = t109 * MDP(20);
t150 = MDP(22) * t159;
t125 = t131 ^ 2;
t149 = t125 * MDP(21) + MDP(18) + 0.2e1 * t150;
t148 = -pkin(5) * t102 - pkin(10) * t101;
t106 = -pkin(5) - t108;
t107 = pkin(10) - t109;
t147 = -t101 * t107 + t102 * t106;
t145 = t121 * MDP(13) + t155;
t144 = MDP(23) * t135 - MDP(24) * t131;
t142 = -MDP(26) * t131 - MDP(27) * t135;
t130 = cos(pkin(6));
t111 = t130 * t137 - t133 * t161;
t112 = t130 * t133 + t137 * t161;
t94 = t129 * t111 - t127 * t112;
t95 = t127 * t111 + t129 * t112;
t88 = t132 * t95 - t136 * t94;
t89 = t132 * t94 + t136 * t95;
t141 = -t89 * MDP(20) - t170 * t88;
t126 = t135 ^ 2;
t85 = t132 * t91 + t136 * t92;
t140 = -t84 * MDP(19) - t85 * MDP(20) + ((-t125 + t126) * MDP(22) + MDP(21) * t159 + MDP(16)) * t102 + (-MDP(17) + t158) * t101;
t83 = t101 * pkin(5) - t102 * pkin(10) + t105;
t79 = t84 * t131;
t78 = -t131 * t160 + t135 * t89;
t77 = -t131 * t89 - t135 * t160;
t76 = t131 * t83 + t135 * t85;
t75 = -t131 * t85 + t135 * t83;
t1 = [MDP(1) + (t128 ^ 2 * t138 ^ 2 + t94 ^ 2 + t95 ^ 2) * MDP(13); (-t94 * t113 - t146 * t95) * MDP(12) + (t94 * t103 + t95 * t104) * MDP(13) + (t77 * t101 + t163 * t88) * MDP(26) + (-t78 * t101 + t162 * t88) * MDP(27) + (-t134 * MDP(4) + (-MDP(11) * t133 - t101 * MDP(19) + MDP(3) - t145 + t157) * t138) * t128; MDP(2) + 0.2e1 * pkin(2) * t157 + (t103 ^ 2 + t104 ^ 2 + t121 ^ 2) * MDP(13) + t155 * t169 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t133 + 0.2e1 * t137 * MDP(6)) * t133 + (t126 * MDP(21) + MDP(14) - 0.2e1 * t150) * t102 ^ 2 + (MDP(19) * t169 + t156 + 0.2e1 * (-MDP(15) + t144) * t102) * t101 + 0.2e1 * (-t103 * t113 - t104 * t146) * MDP(12) + 0.2e1 * (t75 * t101 + t84 * t163) * MDP(26) + 0.2e1 * (-t76 * t101 + t84 * t162) * MDP(27); t111 * MDP(10) - t112 * MDP(11) + (t127 * t95 + t129 * t94) * t165 + t141; (t103 * t129 + t104 * t127) * t165 + t133 * MDP(7) + t137 * MDP(8) + t140 + (t131 * t147 - t164) * MDP(26) + (t135 * t147 + t79) * MDP(27) + (-t129 * t113 - t127 * t146) * pkin(3) * MDP(12) + (-MDP(10) * t133 - MDP(11) * t137) * pkin(8); MDP(9) + (t127 ^ 2 + t129 ^ 2) * MDP(13) * pkin(3) ^ 2 - 0.2e1 * t143 * t106 + 0.2e1 * t154 + 0.2e1 * t153 + t149; -MDP(13) * t160; t170 * t101 + t145; 0; MDP(13); t141; (t131 * t148 - t164) * MDP(26) + (t135 * t148 + t79) * MDP(27) + t140; t149 + t153 + t154 + t143 * (pkin(5) - t106); 0; 0.2e1 * pkin(5) * t143 + t149; t77 * MDP(26) - t78 * MDP(27); t75 * MDP(26) - t76 * MDP(27) + t102 * t144 + t156; t107 * t142 + t158; t143; pkin(10) * t142 + t158; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
