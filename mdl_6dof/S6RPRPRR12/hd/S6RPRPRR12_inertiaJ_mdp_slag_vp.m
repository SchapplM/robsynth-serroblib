% Calculate joint inertia matrix for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR12_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:59
% EndTime: 2019-03-09 04:20:00
% DurationCPUTime: 0.38s
% Computational Cost: add. (390->137), mult. (664->178), div. (0->0), fcn. (602->6), ass. (0->71)
t114 = sin(qJ(6));
t116 = sin(qJ(3));
t117 = cos(qJ(6));
t118 = cos(qJ(5));
t141 = t117 * t118;
t115 = sin(qJ(5));
t143 = t115 * t116;
t83 = t114 * t143 - t116 * t141;
t94 = t114 * t118 + t115 * t117;
t85 = t94 * t116;
t158 = t85 * MDP(27) - t83 * MDP(28);
t120 = -pkin(3) - pkin(8);
t151 = -pkin(9) + t120;
t100 = t151 * t118;
t128 = t114 * t115 - t141;
t98 = t151 * t115;
t130 = -t128 * MDP(27) - t94 * MDP(28) + (t100 * t117 - t114 * t98) * MDP(30) - (t100 * t114 + t117 * t98) * MDP(31);
t157 = (MDP(23) * t120 + MDP(20)) * t118 + (-MDP(24) * t120 - MDP(21)) * t115 + t130;
t126 = MDP(23) * t118 - MDP(24) * t115;
t156 = -2 * MDP(15);
t155 = -2 * MDP(26);
t154 = 0.2e1 * MDP(31);
t153 = (pkin(1) * MDP(6));
t119 = cos(qJ(3));
t152 = t119 * pkin(5);
t84 = t128 * t119;
t86 = t94 * t119;
t150 = t84 * MDP(30) + t86 * MDP(31);
t149 = -MDP(30) * t128 - t94 * MDP(31);
t148 = (MDP(17) * pkin(3));
t97 = t116 * pkin(3) - qJ(4) * t119 + qJ(2);
t88 = pkin(8) * t116 + t97;
t131 = pkin(9) * t116 + t88;
t121 = -pkin(1) - pkin(7);
t101 = (pkin(4) - t121) * t119;
t144 = t115 * t101;
t72 = t131 * t118 + t144;
t147 = t117 * t72;
t146 = MDP(25) * t128;
t145 = MDP(30) * t94;
t142 = t115 * t118;
t110 = t116 ^ 2;
t112 = t119 ^ 2;
t104 = t110 + t112;
t140 = MDP(17) * t121 ^ 2;
t139 = MDP(23) * t115;
t136 = MDP(24) * t118;
t135 = -MDP(13) + MDP(16);
t134 = MDP(22) + MDP(29);
t133 = t119 * MDP(29) + t158;
t132 = MDP(19) * t142;
t93 = t118 * t101;
t71 = -t131 * t115 + t152 + t93;
t68 = -t114 * t72 + t117 * t71;
t129 = MDP(15) - t148;
t127 = MDP(20) * t115 + MDP(21) * t118;
t125 = t136 + t139;
t124 = (MDP(30) * t117 - MDP(31) * t114) * pkin(5);
t123 = MDP(17) * t121 - t126;
t111 = t118 ^ 2;
t109 = t115 ^ 2;
t106 = t116 * t121;
t105 = pkin(5) * t115 + qJ(4);
t102 = pkin(3) * t119 + qJ(4) * t116;
t99 = -pkin(4) * t116 + t106;
t96 = t104 * t121;
t87 = t106 + (-pkin(5) * t118 - pkin(4)) * t116;
t74 = t118 * t88 + t144;
t73 = -t115 * t88 + t93;
t69 = t114 * t71 + t147;
t1 = [MDP(1) + (MDP(17) * t97 + t116 * t156) * t97 + (MDP(25) * t85 + t83 * t155) * t85 + (MDP(18) * t109 + 0.2e1 * t132 + t140) * t110 + ((-2 * MDP(4) + t153) * pkin(1)) + (0.2e1 * MDP(12) * t116 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (MDP(7) + t134 + t140) * t112 + 0.2e1 * (qJ(2) * MDP(13) - t97 * MDP(16) + (-MDP(8) + t127) * t116 + t158) * t119 - 0.2e1 * t96 * MDP(14) + 0.2e1 * (t119 * t68 + t83 * t87) * MDP(30) + (-t119 * t69 + t85 * t87) * t154 + 0.2e1 * (-t119 * t74 + t99 * t143) * MDP(24) + 0.2e1 * (-t116 * t118 * t99 + t119 * t73) * MDP(23); MDP(4) - t153 + t96 * MDP(17) + (t116 * t83 + t119 * t84) * MDP(30) + (t116 * t85 + t119 * t86) * MDP(31) + (-MDP(14) - t126) * t104; MDP(17) * t104 + MDP(6); -t102 * MDP(14) - t85 * t146 + (t128 * t83 - t85 * t94) * MDP(26) + (t105 * t83 + t87 * t94) * MDP(30) + (t105 * t85 - t128 * t87) * MDP(31) + t125 * t99 + (-MDP(10) + MDP(18) * t142 + (-t109 + t111) * MDP(19) + t135 * t121 + t123 * qJ(4)) * t116 + (MDP(9) + (MDP(12) - t129) * t121 + t157) * t119; t102 * MDP(17) + (MDP(12) - MDP(15)) * t119 + (-MDP(31) * t128 + t125 + t135 + t145) * t116; -0.2e1 * t132 + 0.2e1 * t105 * t145 + t111 * MDP(18) + MDP(11) + (t156 + t148) * pkin(3) - (t105 * t154 + t94 * t155 - t146) * t128 + (MDP(17) * qJ(4) + 0.2e1 * MDP(16) + 0.2e1 * t136 + 0.2e1 * t139) * qJ(4); (MDP(14) - t123 + t149) * t119; -t119 * MDP(17); t129; MDP(17); t119 * MDP(22) + t73 * MDP(23) - t74 * MDP(24) + (t117 * t152 + t68) * MDP(30) + (-t147 + (-t71 - t152) * t114) * MDP(31) + t127 * t116 + t133; -t126 * t119 + t150; t157; t126 + t149; 0.2e1 * t124 + t134; MDP(30) * t68 - MDP(31) * t69 + t133; t150; t130; t149; MDP(29) + t124; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
