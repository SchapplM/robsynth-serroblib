% Calculate joint inertia matrix for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP3_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:34
% EndTime: 2019-03-09 08:35:35
% DurationCPUTime: 0.33s
% Computational Cost: add. (364->133), mult. (534->168), div. (0->0), fcn. (403->4), ass. (0->49)
t92 = -pkin(2) - pkin(3);
t88 = sin(qJ(5));
t83 = t88 ^ 2;
t91 = cos(qJ(2));
t115 = t83 * t91;
t89 = sin(qJ(2));
t71 = (pkin(7) - qJ(4)) * t89;
t90 = cos(qJ(5));
t114 = t90 * t71;
t87 = qJ(3) + pkin(4);
t113 = MDP(27) * pkin(5);
t112 = qJ(6) * t91;
t111 = t91 * qJ(3);
t82 = -pkin(8) + t92;
t110 = qJ(6) - t82;
t109 = MDP(22) * t88;
t108 = t88 * MDP(25);
t107 = t90 * MDP(25);
t106 = pkin(7) ^ 2 * MDP(14);
t85 = t90 ^ 2;
t75 = t85 + t83;
t105 = t75 * MDP(27) + MDP(18);
t70 = -t91 * pkin(2) - t89 * qJ(3) - pkin(1);
t103 = t88 * t90 * MDP(20);
t67 = t91 * pkin(3) - t70;
t65 = t89 * pkin(4) + t91 * pkin(8) + t67;
t61 = t90 * t65 - t88 * t71;
t102 = -pkin(2) * MDP(14) - MDP(11);
t101 = MDP(26) * pkin(5) - MDP(21);
t100 = MDP(24) + t113;
t59 = t89 * pkin(5) + t90 * t112 + t61;
t60 = t114 + (t65 + t112) * t88;
t99 = t59 * t90 + t60 * t88;
t98 = -t82 * t89 - t87 * t91;
t97 = t61 * MDP(24) - (t88 * t65 + t114) * MDP(25);
t96 = t90 * MDP(24) - t108;
t95 = -t88 * MDP(24) - MDP(17) - t107;
t93 = qJ(3) ^ 2;
t86 = t91 ^ 2;
t84 = t89 ^ 2;
t79 = t91 * qJ(4);
t77 = t85 * t91;
t76 = t90 * pkin(5) + t87;
t72 = t91 * pkin(7) - t79;
t69 = t110 * t90;
t68 = t110 * t88;
t66 = -t79 + (-pkin(5) * t88 + pkin(7)) * t91;
t63 = -t68 * t88 - t69 * t90;
t1 = [MDP(1) + t70 ^ 2 * MDP(14) + (t67 ^ 2 + t71 ^ 2 + t72 ^ 2) * MDP(18) + (t59 ^ 2 + t60 ^ 2 + t66 ^ 2) * MDP(27) + (t85 * MDP(19) - 0.2e1 * t103 + t106) * t86 + (MDP(23) + MDP(4) + t106) * t84 + 0.2e1 * (-pkin(1) * MDP(10) - t70 * MDP(13) + t67 * MDP(15)) * t89 + 0.2e1 * (-t71 * MDP(17) + t97) * t89 + 0.2e1 * (t84 + t86) * MDP(12) * pkin(7) + 0.2e1 * (-t70 * MDP(11) - t67 * MDP(16) + pkin(1) * MDP(9) + (-MDP(21) * t90 + MDP(5) + t109) * t89 + t99 * MDP(26) + t95 * t72) * t91; t89 * MDP(6) + t91 * MDP(7) + (-t89 * pkin(2) + t111) * MDP(12) + t72 * MDP(15) + t71 * MDP(16) + (-t92 * t89 - t111) * MDP(17) + (t72 * qJ(3) + t71 * t92) * MDP(18) + (t77 - t115) * MDP(20) + (t59 * t68 - t60 * t69 + t66 * t76) * MDP(27) + (-t89 * MDP(22) + t72 * MDP(24) + t98 * MDP(25) + (t68 * t91 - t60) * MDP(26)) * t90 + (t90 * t91 * MDP(19) - t89 * MDP(21) + t98 * MDP(24) - t72 * MDP(25) + (-t69 * t91 + t59) * MDP(26)) * t88 + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t91 + (-MDP(9) + t102) * t89) * pkin(7); MDP(8) + 0.2e1 * pkin(2) * MDP(11) + (pkin(2) ^ 2 + t93) * MDP(14) + (t92 ^ 2 + t93) * MDP(18) + t83 * MDP(19) + 0.2e1 * t103 + (t68 ^ 2 + t69 ^ 2 + t76 ^ 2) * MDP(27) + 0.2e1 * t92 * MDP(16) - 0.2e1 * t63 * MDP(26) + 0.2e1 * t96 * t87 + 0.2e1 * (MDP(13) + MDP(15)) * qJ(3); t71 * MDP(18) + (-t59 * t88 + t60 * t90) * MDP(27) + (pkin(7) * MDP(14) + MDP(12) + t95) * t89; t92 * MDP(18) - t75 * MDP(26) + t63 * MDP(27) + MDP(16) + t102; MDP(14) + t105; -t91 * MDP(16) + t67 * MDP(18) + (t77 + t115) * MDP(26) + t99 * MDP(27) + (MDP(15) + t96) * t89; (t68 * t90 - t69 * t88) * MDP(27); 0; t105; t59 * t113 + t89 * MDP(23) + (t101 * t90 + t109) * t91 + t97; t68 * t113 + (-t82 * MDP(25) - MDP(22)) * t90 + (-t82 * MDP(24) + t101) * t88; -t100 * t88 - t107; t100 * t90 - t108; MDP(27) * pkin(5) ^ 2 + MDP(23); t66 * MDP(27); t76 * MDP(27); 0; 0; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
