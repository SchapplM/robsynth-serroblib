% Calculate joint inertia matrix for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR4_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:48
% EndTime: 2019-03-09 02:26:50
% DurationCPUTime: 0.42s
% Computational Cost: add. (437->127), mult. (743->182), div. (0->0), fcn. (725->8), ass. (0->64)
t124 = sin(qJ(6));
t125 = sin(qJ(5));
t127 = cos(qJ(6));
t128 = cos(qJ(5));
t105 = t124 * t128 + t127 * t125;
t126 = sin(qJ(4));
t95 = t105 * t126;
t161 = -t124 * t125 + t127 * t128;
t96 = t161 * t126;
t164 = -t96 * MDP(26) + t95 * MDP(27);
t134 = -MDP(22) * t125 - MDP(23) * t128;
t158 = pkin(8) + pkin(9);
t109 = t158 * t125;
t110 = t158 * t128;
t138 = t105 * MDP(26) + t161 * MDP(27) + (-t127 * t109 - t124 * t110) * MDP(29) - (-t124 * t109 + t127 * t110) * MDP(30);
t163 = t125 * MDP(19) + t128 * MDP(20) + pkin(8) * t134 + t138;
t122 = sin(pkin(10));
t123 = cos(pkin(10));
t129 = cos(qJ(4));
t148 = t125 * t129;
t98 = -t122 * t148 - t128 * t123;
t146 = t128 * t129;
t99 = t122 * t146 - t125 * t123;
t155 = (-t124 * t99 + t127 * t98) * MDP(29) - (t124 * t98 + t127 * t99) * MDP(30);
t162 = t98 * MDP(22) - t99 * MDP(23) + t155;
t160 = -2 * MDP(25);
t159 = 0.2e1 * MDP(30);
t157 = pkin(9) * t126;
t156 = t129 * pkin(5);
t154 = -t95 * MDP(29) - t96 * MDP(30);
t130 = -pkin(1) - pkin(2);
t108 = t123 * qJ(2) + t122 * t130;
t103 = -pkin(7) + t108;
t140 = t103 * t146;
t106 = t122 * qJ(2) - t123 * t130;
t102 = pkin(3) + t106;
t90 = t129 * pkin(4) + t126 * pkin(8) + t102;
t77 = t140 + (t90 + t157) * t125;
t153 = t127 * t77;
t152 = t103 * t125;
t151 = t103 * t128;
t149 = t125 * t128;
t145 = t161 * MDP(29);
t144 = t105 * MDP(24);
t143 = t126 * MDP(16);
t142 = MDP(21) + MDP(28);
t141 = t129 * MDP(28) + t164;
t139 = MDP(18) * t149;
t88 = t128 * t90;
t76 = t128 * t157 + t88 + (pkin(5) - t152) * t129;
t73 = -t124 * t77 + t127 * t76;
t136 = -t128 * MDP(19) + t125 * MDP(20);
t135 = t128 * MDP(22) - t125 * MDP(23);
t133 = (MDP(29) * t127 - MDP(30) * t124) * pkin(5);
t132 = -t105 * MDP(30) + MDP(15) + t135 + t145;
t120 = t128 ^ 2;
t119 = t126 ^ 2;
t118 = t125 ^ 2;
t115 = -t128 * pkin(5) - pkin(4);
t89 = (-pkin(5) * t125 + t103) * t126;
t79 = t125 * t90 + t140;
t78 = -t103 * t148 + t88;
t74 = t124 * t76 + t153;
t1 = [(t106 ^ 2 + t108 ^ 2) * MDP(9) - 0.2e1 * t102 * t143 + (pkin(1) ^ 2 + qJ(2) ^ 2) * MDP(6) + 0.2e1 * pkin(1) * MDP(4) + 0.2e1 * qJ(2) * MDP(5) + MDP(1) + (t96 * MDP(24) + t95 * t160) * t96 + t142 * t129 ^ 2 + (t120 * MDP(17) + MDP(10) - 0.2e1 * t139) * t119 + 0.2e1 * (t102 * MDP(15) + (MDP(11) + t136) * t126 + t164) * t129 + 0.2e1 * (-t119 * t151 - t79 * t129) * MDP(23) + 0.2e1 * (-t119 * t152 + t78 * t129) * MDP(22) + (-t74 * t129 - t89 * t96) * t159 + 0.2e1 * (t73 * t129 - t89 * t95) * MDP(29) + 0.2e1 * t106 * MDP(7) + 0.2e1 * t108 * MDP(8); -pkin(1) * MDP(6) - MDP(4) + t162 * t129 + (t108 * MDP(9) + t134 * t119 + t154 * t126 + MDP(8)) * t122 + (-t129 * MDP(15) - t106 * MDP(9) - MDP(7) + t143) * t123; MDP(6) + (t122 ^ 2 + t123 ^ 2) * MDP(9); 0; 0; MDP(9); -t96 * t144 + (t105 * t95 - t161 * t96) * MDP(25) + (-t115 * t95 - t161 * t89) * MDP(29) + (t89 * t105 - t115 * t96) * MDP(30) + (-t103 * MDP(16) - MDP(13) + t163) * t129 + (-MDP(12) - t103 * MDP(15) - MDP(17) * t149 + (t118 - t120) * MDP(18) + (pkin(4) * t125 - t151) * MDP(22) + (pkin(4) * t128 + t152) * MDP(23)) * t126; (-t129 * MDP(16) - t126 * t132) * t122; t129 * t132 - t143; 0.2e1 * t139 - 0.2e1 * t115 * t145 + t118 * MDP(17) + MDP(14) + 0.2e1 * t135 * pkin(4) + (t115 * t159 - t160 * t161 + t144) * t105; t129 * MDP(21) + t78 * MDP(22) - t79 * MDP(23) + (t127 * t156 + t73) * MDP(29) + (-t153 + (-t76 - t156) * t124) * MDP(30) + t136 * t126 + t141; t162; t126 * t134 + t154; t163; 0.2e1 * t133 + t142; t73 * MDP(29) - t74 * MDP(30) + t141; t155; t154; t138; MDP(28) + t133; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
