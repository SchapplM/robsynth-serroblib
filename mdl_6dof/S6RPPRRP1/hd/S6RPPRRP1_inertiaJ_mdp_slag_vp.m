% Calculate joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:44
% EndTime: 2019-03-09 01:58:45
% DurationCPUTime: 0.31s
% Computational Cost: add. (498->107), mult. (885->155), div. (0->0), fcn. (911->8), ass. (0->60)
t124 = cos(qJ(4));
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t95 = sin(qJ(4));
t77 = t124 * t90 + t95 * t92;
t127 = 0.2e1 * t77;
t93 = cos(pkin(9));
t84 = -pkin(1) * t93 - pkin(2);
t78 = -pkin(3) * t92 + t84;
t126 = 0.2e1 * t78;
t91 = sin(pkin(9));
t82 = pkin(1) * t91 + qJ(3);
t125 = pkin(7) + t82;
t72 = t125 * t90;
t73 = t125 * t92;
t69 = t124 * t73 - t95 * t72;
t96 = cos(qJ(5));
t123 = t69 * t96;
t94 = sin(qJ(5));
t88 = t94 ^ 2;
t122 = t77 * t88;
t121 = -qJ(6) - pkin(8);
t120 = t90 ^ 2 + t92 ^ 2;
t89 = t96 ^ 2;
t119 = t88 + t89;
t118 = MDP(24) * pkin(5);
t117 = qJ(6) * t77;
t116 = t84 * MDP(8);
t115 = t90 * MDP(6);
t114 = t92 * MDP(5);
t113 = MDP(19) * t94;
t112 = MDP(22) * t94;
t111 = MDP(22) * t96;
t71 = t89 * t77;
t110 = (-t71 - t122) * MDP(23);
t109 = t77 * MDP(15);
t85 = -pkin(5) * t96 - pkin(4);
t108 = t85 * MDP(24);
t107 = t94 * t96 * MDP(17);
t76 = -t124 * t92 + t90 * t95;
t67 = pkin(4) * t76 - pkin(8) * t77 + t78;
t63 = t96 * t67 - t69 * t94;
t106 = t120 * MDP(8);
t105 = -MDP(23) * pkin(5) + MDP(18);
t104 = MDP(21) + t118;
t103 = -pkin(4) * t77 - pkin(8) * t76;
t61 = pkin(5) * t76 - t96 * t117 + t63;
t62 = t123 + (t67 - t117) * t94;
t102 = t61 * t96 + t62 * t94;
t79 = t121 * t94;
t80 = t121 * t96;
t101 = -t79 * t94 - t80 * t96;
t100 = MDP(21) * t63 - MDP(22) * (t67 * t94 + t123);
t99 = MDP(21) * t96 - t112;
t68 = t124 * t72 + t95 * t73;
t98 = MDP(14) + t99;
t75 = t77 ^ 2;
t74 = t76 ^ 2;
t65 = t94 * t77 * pkin(5) + t68;
t1 = [MDP(1) + t109 * t126 + t74 * MDP(20) + (t61 ^ 2 + t62 ^ 2 + t65 ^ 2) * MDP(24) + (t91 ^ 2 + t93 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t114 + 0.2e1 * t115 + t116) * t84 + (MDP(16) * t89 + MDP(9) - 0.2e1 * t107) * t75 + (-t102 * MDP(23) + (MDP(21) * t94 + t111) * t68) * t127 + (0.2e1 * t120 * MDP(7) + t106 * t82) * t82 + (MDP(14) * t126 + (MDP(18) * t96 - MDP(10) - t113) * t127 + 0.2e1 * t100) * t76; (t65 * t76 + (-t61 * t94 + t62 * t96) * t77) * MDP(24); MDP(4) + t106 + (t119 * t75 + t74) * MDP(24); t102 * MDP(24) + t98 * t76 + t109 + t110 - t114 + t115 + t116; 0; t119 * MDP(24) + MDP(8); t77 * MDP(11) - t76 * MDP(12) - t68 * MDP(14) - t69 * MDP(15) + (t71 - t122) * MDP(17) + (t61 * t79 - t62 * t80 + t65 * t85) * MDP(24) + (t76 * MDP(19) - t68 * MDP(21) + t103 * MDP(22) + (-t77 * t79 + t62) * MDP(23)) * t96 + (t96 * t77 * MDP(16) + t76 * MDP(18) + t103 * MDP(21) + t68 * MDP(22) + (t77 * t80 - t61) * MDP(23)) * t94; -t110 + (t101 * MDP(24) - MDP(15)) * t77 + (-t98 + t108) * t76; (t79 * t96 - t80 * t94) * MDP(24); MDP(13) + t88 * MDP(16) + 0.2e1 * t107 + 0.2e1 * t101 * MDP(23) + (t79 ^ 2 + t80 ^ 2 + t85 ^ 2) * MDP(24) + 0.2e1 * t99 * pkin(4); t61 * t118 + MDP(20) * t76 + (t105 * t96 - t113) * t77 + t100; (-t104 * t94 - t111) * t77; t104 * t96 - t112; t79 * t118 + (-MDP(22) * pkin(8) + MDP(19)) * t96 + (-MDP(21) * pkin(8) + t105) * t94; MDP(24) * pkin(5) ^ 2 + MDP(20); t65 * MDP(24); t76 * MDP(24); 0; t108; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
