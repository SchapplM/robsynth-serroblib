% Calculate joint inertia matrix for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR5_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:37
% EndTime: 2019-03-08 20:43:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (342->107), mult. (662->151), div. (0->0), fcn. (708->10), ass. (0->60)
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t138 = t112 * MDP(24) + t116 * MDP(25);
t153 = t112 * MDP(27) + t116 * MDP(28);
t125 = t116 * MDP(27) - t112 * MDP(28);
t152 = -MDP(20) - t125;
t113 = sin(qJ(5));
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t149 = cos(qJ(5));
t93 = t113 * t114 - t149 * t117;
t150 = t93 ^ 2;
t94 = t113 * t117 + t149 * t114;
t91 = t94 ^ 2;
t148 = (pkin(2) * MDP(7));
t111 = cos(pkin(6));
t110 = sin(pkin(6));
t118 = cos(qJ(2));
t141 = t110 * t118;
t85 = -t111 * t114 - t117 * t141;
t86 = t111 * t117 - t114 * t141;
t72 = t113 * t86 - t149 * t85;
t147 = t72 * t93;
t119 = -pkin(2) - pkin(8);
t145 = -pkin(9) + t119;
t96 = t145 * t114;
t97 = t145 * t117;
t79 = t113 * t96 - t149 * t97;
t76 = t79 * t112;
t144 = t79 * t116;
t143 = qJ(3) * MDP(7);
t115 = sin(qJ(2));
t142 = t110 * t115;
t140 = t112 * t116;
t139 = t94 * MDP(20);
t99 = t114 * pkin(4) + qJ(3);
t135 = t114 * MDP(13);
t132 = t117 * MDP(14);
t108 = t112 ^ 2;
t130 = MDP(23) * t140;
t131 = t108 * MDP(22) + MDP(19) + 0.2e1 * t130;
t129 = MDP(5) - t148;
t128 = pkin(5) * t93 - pkin(10) * t94;
t101 = t113 * pkin(4) + pkin(10);
t102 = -t149 * pkin(4) - pkin(5);
t127 = -t101 * t94 - t102 * t93;
t126 = -MDP(24) * t116 + MDP(25) * t112;
t73 = t113 * t85 + t149 * t86;
t123 = -t73 * MDP(21) + t152 * t72;
t122 = -t94 * MDP(21) + t152 * t93;
t109 = t116 ^ 2;
t80 = t113 * t97 + t149 * t96;
t121 = -t79 * MDP(20) - t80 * MDP(21) + (-MDP(18) + t138) * t94 + ((t108 - t109) * MDP(23) - MDP(22) * t140 - MDP(17)) * t93;
t120 = (t149 * MDP(20) - t113 * MDP(21)) * pkin(4);
t75 = t94 * pkin(5) + t93 * pkin(10) + t99;
t68 = t112 * t142 + t116 * t73;
t67 = -t112 * t73 + t116 * t142;
t66 = t112 * t75 + t116 * t80;
t65 = -t112 * t80 + t116 * t75;
t1 = [MDP(1) + (t111 ^ 2 + (t115 ^ 2 + t118 ^ 2) * t110 ^ 2) * MDP(7); (-t112 * t147 + t67 * t94) * MDP(27) + (-t116 * t147 - t68 * t94) * MDP(28) + ((MDP(3) - t129) * t118 + (-t93 * MDP(21) - MDP(4) + MDP(6) + t132 + t135 + t139 + t143) * t115) * t110; 0.2e1 * t99 * t139 + t91 * MDP(26) + MDP(2) + (MDP(8) * t117 - 0.2e1 * t114 * MDP(9)) * t117 + ((-2 * MDP(5) + t148) * pkin(2)) + 0.2e1 * (-t99 * MDP(21) + (MDP(16) + t126) * t94) * t93 + 0.2e1 * (t65 * t94 - t93 * t76) * MDP(27) + 0.2e1 * (-t93 * t144 - t66 * t94) * MDP(28) + (t109 * MDP(22) + MDP(15) - 0.2e1 * t130) * t150 + (0.2e1 * MDP(6) + 0.2e1 * t132 + 0.2e1 * t135 + t143) * qJ(3); -MDP(7) * t141; t129 + t153 * (-t91 - t150); MDP(7); t85 * MDP(13) - t86 * MDP(14) + t123; (t127 * t112 - t144) * MDP(27) + (t127 * t116 + t76) * MDP(28) + (t119 * MDP(13) + MDP(10)) * t117 + (-t119 * MDP(14) - MDP(11)) * t114 + t121; t117 * MDP(13) - t114 * MDP(14) + t122; -0.2e1 * t102 * t125 + MDP(12) + 0.2e1 * t120 + t131; t123; (t128 * t112 - t144) * MDP(27) + (t128 * t116 + t76) * MDP(28) + t121; t122; t120 + t131 + t125 * (pkin(5) - t102); 0.2e1 * pkin(5) * t125 + t131; t67 * MDP(27) - t68 * MDP(28); t94 * MDP(26) + t65 * MDP(27) - t66 * MDP(28) + t126 * t93; -t153 * t94; -t101 * t153 + t138; -pkin(10) * t153 + t138; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
