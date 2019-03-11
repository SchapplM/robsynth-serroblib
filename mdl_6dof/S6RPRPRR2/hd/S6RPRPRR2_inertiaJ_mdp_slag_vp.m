% Calculate joint inertia matrix for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:59
% EndTime: 2019-03-09 03:39:00
% DurationCPUTime: 0.41s
% Computational Cost: add. (668->120), mult. (1228->189), div. (0->0), fcn. (1359->10), ass. (0->68)
t121 = sin(pkin(11));
t123 = cos(pkin(11));
t127 = sin(qJ(3));
t155 = cos(qJ(3));
t108 = t121 * t155 + t123 * t127;
t125 = sin(qJ(6));
t126 = sin(qJ(5));
t128 = cos(qJ(6));
t129 = cos(qJ(5));
t110 = t125 * t129 + t128 * t126;
t86 = t110 * t108;
t109 = t125 * t126 - t128 * t129;
t87 = t109 * t108;
t161 = -t87 * MDP(23) - t86 * MDP(24);
t122 = sin(pkin(10));
t115 = t122 * pkin(1) + pkin(7);
t160 = qJ(4) + t115;
t114 = t121 * pkin(3) + pkin(8);
t136 = t126 * MDP(19) + t129 * MDP(20);
t153 = pkin(9) + t114;
t102 = t153 * t126;
t103 = t153 * t129;
t138 = t110 * MDP(23) - t109 * MDP(24) + (-t128 * t102 - t125 * t103) * MDP(26) - (-t125 * t102 + t128 * t103) * MDP(27);
t159 = t126 * MDP(16) + t129 * MDP(17) - t136 * t114 + t138;
t137 = t129 * MDP(19) - t126 * MDP(20);
t98 = t109 * MDP(26);
t148 = -t110 * MDP(27) - t98;
t133 = t137 + t148;
t106 = t121 * t127 - t123 * t155;
t158 = t106 ^ 2;
t157 = -2 * MDP(22);
t156 = 0.2e1 * MDP(27);
t154 = t106 * pkin(5);
t152 = -t86 * MDP(26) + t87 * MDP(27);
t151 = MDP(13) * pkin(3);
t139 = t160 * t127;
t97 = t160 * t155;
t92 = -t121 * t139 + t123 * t97;
t149 = t129 * t92;
t124 = cos(pkin(10));
t117 = -t124 * pkin(1) - pkin(2);
t112 = -t155 * pkin(3) + t117;
t85 = t106 * pkin(4) - t108 * pkin(8) + t112;
t76 = t149 + (-pkin(9) * t108 + t85) * t126;
t150 = t128 * t76;
t147 = t108 * t126;
t146 = t108 * t129;
t145 = t126 * t129;
t144 = MDP(21) * t110;
t143 = MDP(18) + MDP(25);
t142 = t106 * MDP(25) + t161;
t116 = -t123 * pkin(3) - pkin(4);
t141 = MDP(15) * t145;
t140 = t155 * MDP(10);
t90 = t121 * t97 + t123 * t139;
t77 = -t126 * t92 + t129 * t85;
t75 = -pkin(9) * t146 + t154 + t77;
t72 = -t125 * t76 + t128 * t75;
t135 = (MDP(26) * t128 - MDP(27) * t125) * pkin(5);
t134 = (t129 * MDP(16) - t126 * MDP(17)) * t108;
t120 = t129 ^ 2;
t119 = t126 ^ 2;
t111 = -t129 * pkin(5) + t116;
t105 = t108 ^ 2;
t79 = pkin(5) * t147 + t90;
t78 = t126 * t85 + t149;
t73 = t125 * t75 + t150;
t1 = [MDP(1) - 0.2e1 * t117 * t140 + (t112 ^ 2 + t90 ^ 2 + t92 ^ 2) * MDP(13) - (-t87 * MDP(21) + t86 * t157) * t87 + (t122 ^ 2 + t124 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t120 * MDP(14) - 0.2e1 * t141) * t105 + t143 * t158 + (0.2e1 * t117 * MDP(11) + MDP(5) * t127 + 0.2e1 * t155 * MDP(6)) * t127 + 0.2e1 * (t134 + t161) * t106 + 0.2e1 * (-t92 * t106 + t90 * t108) * MDP(12) + 0.2e1 * (t77 * t106 + t90 * t147) * MDP(19) + 0.2e1 * (-t78 * t106 + t90 * t146) * MDP(20) + 0.2e1 * (t72 * t106 + t79 * t86) * MDP(26) + (-t73 * t106 - t79 * t87) * t156; (t90 * t106 + t92 * t108) * MDP(13); MDP(4) + (t105 + t158) * MDP(13); t127 * MDP(7) + t155 * MDP(8) - t87 * t144 + (t87 * t109 - t110 * t86) * MDP(22) + (t79 * t109 + t111 * t86) * MDP(26) + (t79 * t110 - t111 * t87) * MDP(27) - t137 * t90 + (-t127 * MDP(10) - t155 * MDP(11)) * t115 + (MDP(14) * t145 + (-t119 + t120) * MDP(15) + t136 * t116) * t108 + t159 * t106 + ((-t106 * t121 - t108 * t123) * MDP(12) + (t121 * t92 - t123 * t90) * MDP(13)) * pkin(3); t108 * t121 * t151 + t140 - t127 * MDP(11) + (-t123 * t151 - t133) * t106; 0.2e1 * t141 + 0.2e1 * t111 * t98 + t119 * MDP(14) + MDP(9) + (t121 ^ 2 + t123 ^ 2) * MDP(13) * pkin(3) ^ 2 - 0.2e1 * t137 * t116 + (t109 * t157 + t111 * t156 + t144) * t110; t112 * MDP(13) + t133 * t106; 0; 0; MDP(13); t106 * MDP(18) + t77 * MDP(19) - t78 * MDP(20) + (t128 * t154 + t72) * MDP(26) + (-t150 + (-t75 - t154) * t125) * MDP(27) + t134 + t142; -t136 * t108 + t152; t159; t133; 0.2e1 * t135 + t143; t72 * MDP(26) - t73 * MDP(27) + t142; t152; t138; t148; MDP(25) + t135; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
