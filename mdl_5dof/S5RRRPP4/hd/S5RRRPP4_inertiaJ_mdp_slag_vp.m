% Calculate joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP4_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (503->100), mult. (895->144), div. (0->0), fcn. (922->6), ass. (0->42)
t89 = cos(qJ(2));
t82 = -t89 * pkin(2) - pkin(1);
t107 = 0.2e1 * t82;
t106 = 0.2e1 * t89;
t105 = 2 * MDP(22);
t104 = pkin(6) + pkin(7);
t86 = sin(qJ(3));
t103 = pkin(2) * t86;
t88 = cos(qJ(3));
t81 = t88 * pkin(2) + pkin(3);
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t68 = t85 * t103 + t84 * t81;
t87 = sin(qJ(2));
t72 = t86 * t87 - t88 * t89;
t73 = t86 * t89 + t88 * t87;
t60 = t85 * t72 + t84 * t73;
t61 = -t84 * t72 + t85 * t73;
t64 = t72 * pkin(3) + t82;
t53 = t60 * pkin(4) - t61 * qJ(5) + t64;
t102 = t53 * MDP(23);
t101 = t60 * MDP(20);
t100 = t61 * MDP(22);
t99 = t88 * MDP(16);
t96 = t104 * t89;
t97 = t104 * t87;
t92 = t86 * t97 - t88 * t96;
t59 = -t72 * qJ(4) - t92;
t93 = -t86 * t96 - t88 * t97;
t91 = -t73 * qJ(4) + t93;
t55 = t84 * t59 - t85 * t91;
t57 = t85 * t59 + t84 * t91;
t98 = t55 ^ 2 + t57 ^ 2;
t75 = t85 * t81;
t67 = -t84 * t103 + t75;
t94 = t73 * MDP(13) - t72 * MDP(14) + t93 * MDP(16) + t92 * MDP(17) - t55 * MDP(20) + t57 * MDP(22);
t83 = t84 * pkin(3);
t78 = t85 * pkin(3) + pkin(4);
t77 = t83 + qJ(5);
t66 = -pkin(4) - t67;
t65 = qJ(5) + t68;
t1 = [MDP(1) + pkin(1) * MDP(9) * t106 + t72 * MDP(16) * t107 + (t64 ^ 2 + t98) * MDP(19) + t98 * MDP(23) + (-0.2e1 * t100 + 0.2e1 * t101 + t102) * t53 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t87 + MDP(5) * t106) * t87 + (MDP(11) * t73 - 0.2e1 * t72 * MDP(12) + MDP(17) * t107) * t73 + 0.2e1 * (MDP(18) + MDP(21)) * (t55 * t61 - t57 * t60); t87 * MDP(6) + t89 * MDP(7) + (-t68 * t60 - t67 * t61) * MDP(18) + (-t55 * t67 + t57 * t68) * MDP(19) + (-t65 * t60 + t66 * t61) * MDP(21) + (t55 * t66 + t57 * t65) * MDP(23) + (-t89 * MDP(10) - t87 * MDP(9)) * pkin(6) + t94; MDP(8) + MDP(15) + (t67 ^ 2 + t68 ^ 2) * MDP(19) + (t65 ^ 2 + t66 ^ 2) * MDP(23) + 0.2e1 * (-t86 * MDP(17) + t99) * pkin(2) - 0.2e1 * t66 * MDP(20) + t65 * t105; (-t77 * t60 - t78 * t61) * MDP(21) + (-t55 * t78 + t57 * t77) * MDP(23) + ((-t60 * t84 - t61 * t85) * MDP(18) + (-t55 * t85 + t57 * t84) * MDP(19)) * pkin(3) + t94; MDP(15) + (0.2e1 * pkin(4) + t75) * MDP(20) + (t83 + 0.2e1 * qJ(5) + t68) * MDP(22) + (t65 * t77 - t66 * t78) * MDP(23) + ((t67 * t85 + t68 * t84) * MDP(19) + t85 * MDP(20)) * pkin(3) + (t99 + (-MDP(20) * t84 - MDP(17)) * t86) * pkin(2); MDP(15) + (t77 ^ 2 + t78 ^ 2) * MDP(23) + (t84 ^ 2 + t85 ^ 2) * MDP(19) * pkin(3) ^ 2 + 0.2e1 * t78 * MDP(20) + t77 * t105; t64 * MDP(19) - t100 + t101 + t102; 0; 0; MDP(19) + MDP(23); t61 * MDP(21) + t55 * MDP(23); t66 * MDP(23) - MDP(20); -t78 * MDP(23) - MDP(20); 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
