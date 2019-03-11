% Calculate joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR2_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:23
% EndTime: 2019-03-09 02:21:24
% DurationCPUTime: 0.36s
% Computational Cost: add. (609->112), mult. (1140->163), div. (0->0), fcn. (1270->10), ass. (0->71)
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t133 = -t124 * MDP(21) - t127 * MDP(22);
t123 = sin(qJ(6));
t126 = cos(qJ(6));
t104 = t123 * t124 - t126 * t127;
t105 = t123 * t127 + t126 * t124;
t157 = pkin(8) + pkin(9);
t107 = t157 * t124;
t108 = t157 * t127;
t137 = t105 * MDP(25) - t104 * MDP(26) + (-t126 * t107 - t123 * t108) * MDP(28) - (-t123 * t107 + t126 * t108) * MDP(29);
t162 = t124 * MDP(18) + t127 * MDP(19) + t133 * pkin(8) + t137;
t119 = sin(pkin(11));
t121 = cos(pkin(11));
t125 = sin(qJ(4));
t156 = cos(qJ(4));
t102 = t125 * t119 - t121 * t156;
t134 = t127 * MDP(21) - t124 * MDP(22);
t132 = -MDP(14) - t134;
t103 = t119 * t156 + t125 * t121;
t142 = t103 * MDP(15);
t96 = t104 * MDP(28);
t152 = -t105 * MDP(29) - t96;
t161 = -(-t132 + t152) * t102 - t142;
t122 = cos(pkin(10));
t112 = -t122 * pkin(1) - pkin(2);
t106 = -t121 * pkin(3) + t112;
t160 = 0.2e1 * t106;
t159 = -2 * MDP(24);
t158 = 0.2e1 * MDP(29);
t155 = t102 * pkin(5);
t120 = sin(pkin(10));
t110 = t120 * pkin(1) + qJ(3);
t154 = pkin(7) + t110;
t85 = t105 * t103;
t86 = t104 * t103;
t153 = -t85 * MDP(28) + t86 * MDP(29);
t93 = t154 * t119;
t94 = t154 * t121;
t88 = -t125 * t93 + t156 * t94;
t150 = t127 * t88;
t84 = t102 * pkin(4) - t103 * pkin(8) + t106;
t75 = t150 + (-pkin(9) * t103 + t84) * t124;
t151 = t126 * t75;
t149 = t103 * t124;
t148 = t103 * t127;
t147 = t112 * MDP(8);
t146 = t119 * MDP(6);
t145 = t121 * MDP(5);
t144 = t124 * t127;
t81 = t85 * MDP(26);
t82 = t86 * MDP(25);
t143 = t119 ^ 2 + t121 ^ 2;
t141 = t105 * MDP(23);
t140 = MDP(20) + MDP(27);
t139 = t102 * MDP(27) - t81 - t82;
t138 = MDP(17) * t144;
t76 = -t124 * t88 + t127 * t84;
t74 = -pkin(9) * t148 + t155 + t76;
t71 = -t123 * t75 + t126 * t74;
t136 = t143 * MDP(8);
t135 = t127 * MDP(18) - t124 * MDP(19);
t87 = t125 * t94 + t156 * t93;
t131 = (MDP(28) * t126 - MDP(29) * t123) * pkin(5);
t118 = t127 ^ 2;
t117 = t124 ^ 2;
t114 = -t127 * pkin(5) - pkin(4);
t78 = pkin(5) * t149 + t87;
t77 = t124 * t84 + t150;
t72 = t123 * t74 + t151;
t1 = [t142 * t160 + MDP(1) - (-t86 * MDP(23) + t85 * t159) * t86 + (t120 ^ 2 + t122 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t145 + 0.2e1 * t146 + t147) * t112 + t140 * t102 ^ 2 + (t118 * MDP(16) + MDP(9) - 0.2e1 * t138) * t103 ^ 2 + (MDP(14) * t160 - 0.2e1 * t82 - 0.2e1 * t81 + 0.2e1 * (-MDP(10) + t135) * t103) * t102 + 0.2e1 * (t76 * t102 + t87 * t149) * MDP(21) + 0.2e1 * (-t77 * t102 + t87 * t148) * MDP(22) + 0.2e1 * (t71 * t102 + t78 * t85) * MDP(28) + (-t72 * t102 - t78 * t86) * t158 + (0.2e1 * t143 * MDP(7) + t136 * t110) * t110; 0; MDP(4) + t136; -t145 + t146 + t147 - t161; 0; MDP(8); -t88 * MDP(15) - t86 * t141 + (t86 * t104 - t105 * t85) * MDP(24) + (t78 * t104 + t114 * t85) * MDP(28) + (t78 * t105 - t114 * t86) * MDP(29) + t132 * t87 + (MDP(11) + MDP(16) * t144 + (-t117 + t118) * MDP(17) + t133 * pkin(4)) * t103 + (-MDP(12) + t162) * t102; t161; 0; 0.2e1 * t138 + 0.2e1 * t114 * t96 + t117 * MDP(16) + MDP(13) + 0.2e1 * t134 * pkin(4) + (t104 * t159 + t114 * t158 + t141) * t105; t102 * MDP(20) + t76 * MDP(21) - t77 * MDP(22) + (t126 * t155 + t71) * MDP(28) + (-t151 + (-t74 - t155) * t123) * MDP(29) + t135 * t103 + t139; t133 * t103 + t153; t134 + t152; t162; 0.2e1 * t131 + t140; t71 * MDP(28) - t72 * MDP(29) + t139; t153; t152; t137; MDP(27) + t131; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
