% Calculate joint inertia matrix for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP13_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:38
% EndTime: 2019-12-31 18:59:39
% DurationCPUTime: 0.31s
% Computational Cost: add. (261->105), mult. (455->140), div. (0->0), fcn. (339->4), ass. (0->41)
t89 = pkin(7) * MDP(24);
t98 = MDP(22) + t89;
t67 = sin(qJ(4));
t69 = cos(qJ(4));
t85 = -MDP(20) + MDP(23);
t86 = MDP(19) + MDP(21);
t97 = -t86 * t67 + t85 * t69;
t96 = 2 * pkin(4);
t70 = cos(qJ(3));
t95 = 0.2e1 * t70;
t94 = (pkin(1) * MDP(6));
t71 = -pkin(1) - pkin(6);
t93 = t67 * t71;
t92 = t69 * t71;
t68 = sin(qJ(3));
t60 = t68 * pkin(3) - t70 * pkin(7) + qJ(2);
t54 = t67 * t60 + t68 * t92;
t63 = t67 ^ 2;
t65 = t69 ^ 2;
t91 = t63 + t65;
t80 = -t69 * pkin(4) - t67 * qJ(5);
t61 = -pkin(3) + t80;
t88 = t61 * MDP(24);
t87 = t69 * MDP(15);
t84 = t91 * MDP(22);
t83 = -MDP(24) * pkin(4) - MDP(21);
t81 = 0.2e1 * qJ(5) * MDP(23) + MDP(18);
t79 = -pkin(4) * t67 + t69 * qJ(5);
t51 = t68 * qJ(5) + t54;
t57 = t69 * t60;
t52 = -t57 + (-pkin(4) + t93) * t68;
t78 = t51 * t69 + t52 * t67;
t77 = t69 * MDP(16) - t67 * MDP(17);
t76 = t67 * MDP(16) + t69 * MDP(17);
t75 = (-t68 * t93 + t57) * MDP(19) - t54 * MDP(20);
t74 = t67 * MDP(21) - t69 * MDP(23);
t73 = (MDP(24) * qJ(5) + t85) * t69 + (-MDP(19) + t83) * t67;
t66 = t70 ^ 2;
t64 = t68 ^ 2;
t55 = (-t71 - t79) * t70;
t1 = [MDP(1) + t64 * MDP(18) + (t51 ^ 2 + t52 ^ 2 + t55 ^ 2) * MDP(24) + ((-2 * MDP(4) + t94) * pkin(1)) + (MDP(13) * t95 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + ((-t51 * t67 + t52 * t69) * MDP(22) + t74 * t55) * t95 + (t65 * MDP(14) - 0.2e1 * t67 * t87 + MDP(7) + 0.2e1 * (-t67 * MDP(19) - t69 * MDP(20)) * t71) * t66 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t77) * t70 - t52 * MDP(21) + t51 * MDP(23) + t75) * t68; MDP(4) - t94 + (-t55 * t70 + t78 * t68) * MDP(24) + t97 * (t64 + t66); MDP(6) + (t91 * t64 + t66) * MDP(24); (-t69 * MDP(21) - t67 * MDP(23) + t88) * t55 + (-t71 * MDP(13) + t97 * pkin(7) - MDP(10) + t76) * t68 + (MDP(9) + t71 * MDP(12) + t69 * t67 * MDP(14) + (-t63 + t65) * MDP(15) + (-pkin(3) * t67 + t92) * MDP(19) + (-pkin(3) * t69 - t93) * MDP(20) + t74 * t61) * t70 + t98 * t78; (t91 * t89 - MDP(13) + t84) * t68 + (t85 * t67 + t86 * t69 + MDP(12) - t88) * t70; MDP(11) + t63 * MDP(14) + (t91 * pkin(7) ^ 2 + t61 ^ 2) * MDP(24) + 0.2e1 * pkin(7) * t84 + 0.2e1 * (pkin(3) * MDP(19) - t61 * MDP(21)) * t69 + 0.2e1 * (-pkin(3) * MDP(20) - t61 * MDP(23) + t87) * t67; t57 * MDP(21) + t54 * MDP(23) + (-t52 * pkin(4) + t51 * qJ(5)) * MDP(24) + ((t96 - t93) * MDP(21) + t81) * t68 + (t80 * MDP(22) + t77) * t70 + t75; t73 * t68; t79 * MDP(22) + t73 * pkin(7) + t76; MDP(21) * t96 + ((pkin(4) ^ 2) + qJ(5) ^ 2) * MDP(24) + t81; t69 * t70 * MDP(22) - t68 * MDP(21) + t52 * MDP(24); t67 * t68 * MDP(24); t98 * t67; t83; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
