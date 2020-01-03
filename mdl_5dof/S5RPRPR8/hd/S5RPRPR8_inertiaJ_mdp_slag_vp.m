% Calculate joint inertia matrix for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR8_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:12
% EndTime: 2019-12-31 18:22:13
% DurationCPUTime: 0.32s
% Computational Cost: add. (284->91), mult. (553->140), div. (0->0), fcn. (509->8), ass. (0->51)
t117 = qJ(4) * MDP(15);
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t108 = t84 ^ 2 + t86 ^ 2;
t116 = t108 * t117;
t114 = -2 * MDP(17);
t113 = 2 * MDP(22);
t89 = sin(qJ(3));
t112 = pkin(7) * t89;
t85 = sin(pkin(8));
t77 = t85 * pkin(1) + pkin(6);
t111 = t77 * t84;
t91 = cos(qJ(3));
t110 = t77 * t91;
t109 = pkin(7) + qJ(4);
t87 = cos(pkin(8));
t78 = -t87 * pkin(1) - pkin(2);
t70 = -t91 * pkin(3) - t89 * qJ(4) + t78;
t62 = t86 * t110 + t84 * t70;
t107 = pkin(3) * MDP(15);
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t71 = t88 * t84 - t90 * t86;
t66 = t71 * t89;
t106 = t66 * MDP(16);
t105 = t71 * MDP(21);
t104 = t84 * MDP(13);
t103 = t86 * MDP(12);
t102 = t108 * MDP(14);
t68 = t86 * t70;
t61 = -t84 * t110 + t68;
t100 = -t61 * t84 + t62 * t86;
t72 = t90 * t84 + t88 * t86;
t99 = MDP(12) * t84 + MDP(13) * t86;
t65 = t72 * t89;
t98 = t66 * MDP(18) + t65 * MDP(19);
t97 = -t65 * MDP(21) + t66 * MDP(22);
t96 = -t103 + t104 - t107;
t74 = t109 * t84;
t75 = t109 * t86;
t95 = t72 * MDP(18) - t71 * MDP(19) + (-t90 * t74 - t88 * t75) * MDP(21) - (-t88 * t74 + t90 * t75) * MDP(22);
t94 = t72 * MDP(22) + t105 + t96;
t83 = t91 ^ 2;
t82 = t89 ^ 2;
t79 = -t86 * pkin(4) - pkin(3);
t69 = (pkin(4) * t84 + t77) * t89;
t60 = -t84 * t112 + t62;
t59 = -t86 * t112 + t68 + (-pkin(4) - t111) * t91;
t58 = t88 * t59 + t90 * t60;
t57 = t90 * t59 - t88 * t60;
t1 = [MDP(1) + t82 * MDP(5) + (t82 * t77 ^ 2 + t61 ^ 2 + t62 ^ 2) * MDP(15) + t83 * MDP(20) + (t85 ^ 2 + t87 ^ 2) * MDP(4) * pkin(1) ^ 2 - (t65 * t114 - t106) * t66 + 0.2e1 * (-t78 * MDP(10) + t89 * MDP(6) + t98) * t91 + 0.2e1 * (t82 * t111 - t61 * t91) * MDP(12) + 0.2e1 * (t82 * t77 * t86 + t62 * t91) * MDP(13) + 0.2e1 * (-t57 * t91 + t69 * t65) * MDP(21) + (t58 * t91 - t69 * t66) * t113 + 0.2e1 * (t78 * MDP(11) + (-t61 * t86 - t62 * t84) * MDP(14)) * t89; (t100 - t110) * t89 * MDP(15); MDP(4) + (t108 * t82 + t83) * MDP(15); -t72 * t106 + (-t72 * t65 + t66 * t71) * MDP(17) + (t79 * t65 + t69 * t71) * MDP(21) + (-t79 * t66 + t69 * t72) * MDP(22) + (-t77 * MDP(11) + t99 * qJ(4) + MDP(8) - t95) * t91 + (MDP(7) - t99 * pkin(3) + (-MDP(10) + t96) * t77) * t89 + (MDP(14) + t117) * t100; (-MDP(11) + t102 + t116) * t89 + (MDP(10) - t94) * t91; 0.2e1 * t79 * t105 + MDP(9) + (0.2e1 * t103 - 0.2e1 * t104 + t107) * pkin(3) + (MDP(16) * t72 + t79 * t113 + t71 * t114) * t72 + (0.2e1 * t102 + t116) * qJ(4); (t77 * MDP(15) + t99) * t89 - t97; -t91 * MDP(15); t94; MDP(15); -t91 * MDP(20) + t57 * MDP(21) - t58 * MDP(22) - t98; t97; t95; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
