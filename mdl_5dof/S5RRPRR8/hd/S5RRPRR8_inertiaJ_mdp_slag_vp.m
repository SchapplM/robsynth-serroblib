% Calculate joint inertia matrix for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR8_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:14
% EndTime: 2019-12-31 20:18:15
% DurationCPUTime: 0.29s
% Computational Cost: add. (499->93), mult. (950->137), div. (0->0), fcn. (1063->8), ass. (0->53)
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t125 = t101 * MDP(22) + t104 * MDP(23);
t110 = t104 * MDP(25) - t101 * MDP(26);
t100 = cos(pkin(9));
t103 = sin(qJ(2));
t106 = cos(qJ(2));
t99 = sin(pkin(9));
t112 = t100 * t106 - t99 * t103;
t93 = -t106 * pkin(2) - pkin(1);
t79 = -t112 * pkin(3) + t93;
t130 = 0.2e1 * t79;
t129 = 0.2e1 * t106;
t128 = t99 * pkin(2);
t126 = -qJ(3) - pkin(6);
t89 = t126 * t103;
t90 = t126 * t106;
t78 = -t100 * t90 + t99 * t89;
t102 = sin(qJ(4));
t105 = cos(qJ(4));
t77 = t100 * t89 + t99 * t90;
t85 = t100 * t103 + t99 * t106;
t67 = -t85 * pkin(7) + t77;
t68 = t112 * pkin(7) + t78;
t64 = t102 * t68 - t105 * t67;
t60 = t64 * t101;
t124 = t64 * t104;
t123 = t101 * t104;
t75 = t102 * t85 - t105 * t112;
t122 = t75 * MDP(24);
t76 = t102 * t112 + t105 * t85;
t121 = t76 * MDP(19);
t92 = t100 * pkin(2) + pkin(3);
t82 = -t102 * t128 + t105 * t92;
t120 = t82 * MDP(18);
t83 = -t102 * t92 - t105 * t128;
t119 = t83 * MDP(19);
t115 = MDP(21) * t123;
t97 = t101 ^ 2;
t116 = t97 * MDP(20) + MDP(17) + 0.2e1 * t115;
t114 = -pkin(4) * t76 - pkin(8) * t75;
t80 = -pkin(4) - t82;
t81 = pkin(8) - t83;
t113 = -t75 * t81 + t76 * t80;
t111 = MDP(22) * t104 - MDP(23) * t101;
t109 = -MDP(25) * t101 - MDP(26) * t104;
t65 = t102 * t67 + t105 * t68;
t98 = t104 ^ 2;
t108 = -t64 * MDP(18) - t65 * MDP(19) + ((-t97 + t98) * MDP(21) + MDP(20) * t123 + MDP(15)) * t76 + (-MDP(16) + t125) * t75;
t63 = t75 * pkin(4) - t76 * pkin(8) + t79;
t59 = t101 * t63 + t104 * t65;
t58 = -t101 * t65 + t104 * t63;
t1 = [MDP(1) + pkin(1) * MDP(9) * t129 + (t77 ^ 2 + t78 ^ 2 + t93 ^ 2) * MDP(12) + t121 * t130 + (t98 * MDP(20) + MDP(13) - 0.2e1 * t115) * t76 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t103 + MDP(5) * t129) * t103 + (MDP(18) * t130 + t122 + 0.2e1 * (-MDP(14) + t111) * t76) * t75 + 0.2e1 * (t78 * t112 - t77 * t85) * MDP(11) + 0.2e1 * (t58 * t75 + t76 * t60) * MDP(25) + 0.2e1 * (t76 * t124 - t59 * t75) * MDP(26); (t113 * t104 + t60) * MDP(26) + (t113 * t101 - t124) * MDP(25) + t103 * MDP(6) + t106 * MDP(7) + t108 + (-t106 * MDP(10) - t103 * MDP(9)) * pkin(6) + ((-t100 * t85 + t99 * t112) * MDP(11) + (t100 * t77 + t78 * t99) * MDP(12)) * pkin(2); MDP(8) - 0.2e1 * t110 * t80 + (t100 ^ 2 + t99 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t120 + 0.2e1 * t119 + t116; t93 * MDP(12) + t121 + (MDP(18) + t110) * t75; 0; MDP(12); (t114 * t101 - t124) * MDP(25) + (t114 * t104 + t60) * MDP(26) + t108; t116 + t119 + t120 + t110 * (pkin(4) - t80); 0; 0.2e1 * pkin(4) * t110 + t116; t58 * MDP(25) - t59 * MDP(26) + t111 * t76 + t122; t109 * t81 + t125; t110; t109 * pkin(8) + t125; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
