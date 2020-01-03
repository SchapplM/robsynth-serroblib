% Calculate joint inertia matrix for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP7_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:56
% EndTime: 2019-12-31 21:05:57
% DurationCPUTime: 0.50s
% Computational Cost: add. (387->146), mult. (686->190), div. (0->0), fcn. (534->4), ass. (0->47)
t89 = sin(qJ(2));
t118 = 0.2e1 * t89;
t90 = cos(qJ(3));
t112 = qJ(4) * t90;
t88 = sin(qJ(3));
t92 = pkin(3) + pkin(4);
t117 = t88 * t92 - t112;
t91 = cos(qJ(2));
t116 = t88 * t91;
t115 = pkin(7) - qJ(5);
t76 = -pkin(2) * t91 - pkin(7) * t89 - pkin(1);
t114 = pkin(6) * t116 - t90 * t76;
t70 = t90 * t91 * pkin(6) + t88 * t76;
t85 = t88 ^ 2;
t87 = t90 ^ 2;
t113 = t85 + t87;
t111 = qJ(5) * t89;
t110 = t91 * qJ(4);
t109 = MDP(12) * t90;
t108 = MDP(17) - MDP(20);
t107 = MDP(19) - MDP(24);
t84 = t91 * pkin(3);
t68 = t84 + t114;
t106 = -0.2e1 * t110 + t70;
t105 = qJ(4) * t88 + pkin(2);
t104 = -pkin(3) * MDP(21) - MDP(18);
t67 = -t110 + t70;
t103 = -pkin(3) * t88 + t112;
t102 = t67 * t90 + t68 * t88;
t101 = -MDP(16) * t114 - t70 * MDP(17);
t100 = -t88 * MDP(22) + t90 * MDP(23);
t77 = t115 * t88;
t78 = t115 * t90;
t99 = -t88 * MDP(13) + t77 * MDP(22) - t78 * MDP(23);
t72 = t90 * t92 + t105;
t98 = t90 * MDP(22) + t88 * MDP(23) + MDP(25) * t72;
t75 = -pkin(3) * t90 - t105;
t97 = pkin(2) * MDP(16) - t75 * MDP(18) + t72 * MDP(22) - t78 * MDP(24);
t96 = -MDP(17) * pkin(2) - MDP(20) * t75 + MDP(23) * t72 - t77 * MDP(24);
t94 = qJ(4) ^ 2;
t81 = pkin(7) * t116;
t79 = t88 * t111;
t71 = (pkin(6) - t103) * t89;
t66 = (-pkin(6) - t117) * t89;
t65 = t67 + t79;
t64 = pkin(4) * t91 - t90 * t111 + t68;
t1 = [MDP(1) + (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) * MDP(21) + (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) * MDP(25) + (t91 * MDP(15) + 0.2e1 * pkin(1) * MDP(9) + (-t90 * MDP(13) + t88 * MDP(14) + MDP(5)) * t118) * t91 + 0.2e1 * (t68 * MDP(18) - t67 * MDP(20) + t64 * MDP(22) - t65 * MDP(23) - t101) * t91 + ((-t67 * t88 + t68 * t90) * MDP(19) + (-t64 * t90 + t65 * t88) * MDP(24) + (t88 * MDP(18) - t90 * MDP(20)) * t71 + t100 * t66) * t118 + (-0.2e1 * pkin(1) * MDP(10) + (t87 * MDP(11) - 0.2e1 * t88 * t109 + MDP(4) + 0.2e1 * (t88 * MDP(16) + t90 * MDP(17)) * pkin(6)) * t89) * t89; t81 * MDP(16) + (-t71 * t90 + t81) * MDP(18) + t102 * MDP(19) - t71 * t88 * MDP(20) + (t102 * pkin(7) + t71 * t75) * MDP(21) + (-t64 * t88 - t65 * t90) * MDP(24) + (t64 * t77 + t65 * t78) * MDP(25) + t98 * t66 + (-pkin(6) * MDP(10) + MDP(7) + (t108 * pkin(7) - MDP(14)) * t90 + t99) * t91 + (MDP(6) - pkin(6) * MDP(9) + (-t85 + t87) * MDP(12) + (-pkin(6) * MDP(16) + t96) * t90 + (t90 * MDP(11) + pkin(6) * MDP(17) - t97) * t88) * t89; MDP(8) + t85 * MDP(11) + (t113 * pkin(7) ^ 2 + t75 ^ 2) * MDP(21) + (t72 ^ 2 + t77 ^ 2 + t78 ^ 2) * MDP(25) + 0.2e1 * t113 * MDP(19) * pkin(7) + 0.2e1 * t97 * t90 + 0.2e1 * (t96 + t109) * t88; (-0.2e1 * t84 - t114) * MDP(18) + t106 * MDP(20) + (-pkin(3) * t68 + qJ(4) * t67) * MDP(21) - t68 * MDP(22) + (t79 + t106) * MDP(23) + (qJ(4) * t65 - t64 * t92) * MDP(25) + (-MDP(15) + (-pkin(4) - t92) * MDP(22)) * t91 + ((-t107 * qJ(4) - MDP(14)) * t88 + (-pkin(3) * MDP(19) + qJ(5) * MDP(22) + t92 * MDP(24) + MDP(13)) * t90) * t89 + t101; t90 * MDP(14) + t103 * MDP(19) + t117 * MDP(24) + (qJ(4) * t78 - t77 * t92) * MDP(25) + ((MDP(21) * qJ(4) - t108) * t90 + (-MDP(16) + t104) * t88) * pkin(7) - t99; MDP(15) + 0.2e1 * pkin(3) * MDP(18) + (pkin(3) ^ 2 + t94) * MDP(21) + 0.2e1 * t92 * MDP(22) + (t92 ^ 2 + t94) * MDP(25) + 0.2e1 * (MDP(20) + MDP(23)) * qJ(4); t68 * MDP(21) + t64 * MDP(25) + (MDP(18) + MDP(22)) * t91 + t107 * t90 * t89; t77 * MDP(25) + (MDP(21) * pkin(7) + t107) * t88; -MDP(25) * t92 - MDP(22) + t104; MDP(21) + MDP(25); t66 * MDP(25) + t100 * t89; t98; 0; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
