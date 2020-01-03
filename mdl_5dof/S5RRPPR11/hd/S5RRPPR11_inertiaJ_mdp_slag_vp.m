% Calculate joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR11_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:00
% EndTime: 2019-12-31 19:48:01
% DurationCPUTime: 0.32s
% Computational Cost: add. (312->103), mult. (558->144), div. (0->0), fcn. (491->6), ass. (0->50)
t126 = pkin(3) + pkin(6);
t101 = cos(qJ(5));
t100 = sin(qJ(2));
t102 = cos(qJ(2));
t115 = -t100 * qJ(3) - pkin(1);
t98 = -pkin(2) - qJ(4);
t75 = t102 * t98 + t115;
t84 = t126 * t100;
t97 = cos(pkin(8));
t80 = t97 * t84;
t96 = sin(pkin(8));
t65 = t100 * pkin(4) + t80 + (pkin(7) * t102 - t75) * t96;
t119 = t102 * t97;
t68 = t97 * t75 + t96 * t84;
t66 = -pkin(7) * t119 + t68;
t99 = sin(qJ(5));
t120 = t99 * t96;
t71 = -t101 * t119 + t102 * t120;
t77 = t101 * t96 + t99 * t97;
t72 = t77 * t102;
t125 = (t101 * t65 - t99 * t66) * MDP(24) - (t101 * t66 + t99 * t65) * MDP(25) - t72 * MDP(21) + t71 * MDP(22);
t88 = t96 * pkin(4) + qJ(3);
t122 = 0.2e1 * t88;
t121 = -pkin(7) + t98;
t85 = t126 * t102;
t87 = t96 ^ 2 + t97 ^ 2;
t67 = -t96 * t75 + t80;
t64 = t67 * t97 + t68 * t96;
t118 = t64 * MDP(18);
t117 = t72 * MDP(19);
t116 = t77 * MDP(24);
t114 = -pkin(2) * MDP(14) + MDP(12);
t113 = MDP(15) * t97 - MDP(16) * t96;
t112 = t96 * MDP(15) + t97 * MDP(16);
t109 = -t71 * MDP(24) - t72 * MDP(25);
t78 = t101 * t97 - t120;
t108 = t78 * MDP(24) - t77 * MDP(25);
t107 = MDP(11) + t113;
t106 = qJ(3) * MDP(18) + t112;
t81 = t121 * t96;
t82 = t121 * t97;
t105 = t78 * MDP(21) - t77 * MDP(22) + (t101 * t82 - t99 * t81) * MDP(24) - (t101 * t81 + t99 * t82) * MDP(25);
t104 = pkin(6) ^ 2;
t103 = qJ(3) ^ 2;
t95 = t102 ^ 2;
t94 = t100 ^ 2;
t83 = -t102 * pkin(2) + t115;
t76 = t87 * t98;
t74 = pkin(4) * t119 + t85;
t1 = [MDP(1) + (t95 * t104 + t83 ^ 2) * MDP(14) + (t67 ^ 2 + t68 ^ 2 + t85 ^ 2) * MDP(18) + (t104 * MDP(14) + MDP(23) + MDP(4)) * t94 - (0.2e1 * t71 * MDP(20) - t117) * t72 + 0.2e1 * t109 * t74 + 0.2e1 * (t94 + t95) * MDP(11) * pkin(6) + 0.2e1 * (-pkin(1) * MDP(10) - t83 * MDP(13) + t67 * MDP(15) - t68 * MDP(16) + t102 * MDP(5) + t125) * t100 + 0.2e1 * ((t67 * t96 - t68 * t97) * MDP(17) + t113 * t85 + t83 * MDP(12) + pkin(1) * MDP(9)) * t102; -t64 * MDP(17) - t78 * t117 + (t78 * t71 + t72 * t77) * MDP(20) + (-t88 * t71 + t74 * t77) * MDP(24) + (-t88 * t72 + t74 * t78) * MDP(25) + t98 * t118 + t106 * t85 + (qJ(3) * t107 + MDP(7)) * t102 + (-pkin(2) * MDP(11) + t113 * t98 + MDP(6) + t105) * t100 + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t102 + (-MDP(9) + t114) * t100) * pkin(6); MDP(8) - 0.2e1 * pkin(2) * MDP(12) + (pkin(2) ^ 2 + t103) * MDP(14) - 0.2e1 * t76 * MDP(17) + (t87 * t98 ^ 2 + t103) * MDP(18) + t116 * t122 + (MDP(19) * t78 - 0.2e1 * t77 * MDP(20) + MDP(25) * t122) * t78 + 0.2e1 * (MDP(13) + t112) * qJ(3); t118 + (MDP(14) * pkin(6) + t107 + t108) * t100; -t87 * MDP(17) + t76 * MDP(18) + t114; t87 * MDP(18) + MDP(14); t85 * MDP(18) + t102 * t113 + t109; t78 * MDP(25) + t106 + t116; 0; MDP(18); t100 * MDP(23) + t125; t105; t108; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
