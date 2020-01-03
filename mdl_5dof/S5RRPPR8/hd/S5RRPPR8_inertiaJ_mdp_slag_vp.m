% Calculate joint inertia matrix for
% S5RRPPR8
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
%   see S5RRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR8_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:16
% EndTime: 2019-12-31 19:39:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (269->93), mult. (451->136), div. (0->0), fcn. (423->6), ass. (0->36)
t89 = sin(qJ(2));
t91 = cos(qJ(2));
t77 = -t91 * pkin(2) - t89 * qJ(3) - pkin(1);
t65 = t91 * pkin(3) - t77;
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t67 = t89 * t86 + t91 * t87;
t106 = 0.2e1 * t67 * pkin(4) + 0.2e1 * t65;
t105 = -pkin(2) - pkin(3);
t104 = pkin(6) - qJ(4);
t78 = t104 * t91;
t96 = t104 * t89;
t63 = t87 * t78 + t86 * t96;
t84 = t89 ^ 2;
t103 = t91 ^ 2 + t84;
t68 = -t91 * t86 + t89 * t87;
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t56 = t90 * t67 + t88 * t68;
t102 = t56 * MDP(24);
t74 = t86 * qJ(3) - t87 * t105;
t73 = -pkin(4) - t74;
t76 = t87 * qJ(3) + t86 * t105;
t101 = (-t90 * t73 + t88 * t76) * MDP(24);
t100 = (t88 * t73 + t90 * t76) * MDP(25);
t99 = t65 * MDP(18);
t98 = t67 * MDP(15);
t97 = t68 * MDP(16);
t61 = t86 * t78 - t87 * t96;
t95 = -pkin(2) * MDP(14) - MDP(11);
t94 = -(t88 * t86 - t90 * t87) * MDP(24) - (t90 * t86 + t88 * t87) * MDP(25);
t54 = -t68 * pkin(7) - t61;
t55 = -t67 * pkin(7) + t63;
t57 = -t88 * t67 + t90 * t68;
t93 = t57 * MDP(21) - t56 * MDP(22) - (-t90 * t54 + t88 * t55) * MDP(24) - (t88 * t54 + t90 * t55) * MDP(25);
t1 = [MDP(1) + t84 * MDP(4) + (t103 * pkin(6) ^ 2 + t77 ^ 2) * MDP(14) + (t61 ^ 2 + t63 ^ 2) * MDP(18) + t102 * t106 + (0.2e1 * t97 + 0.2e1 * t98 + t99) * t65 + 0.2e1 * (t61 * t68 - t63 * t67) * MDP(17) + 0.2e1 * t103 * MDP(12) * pkin(6) + (MDP(19) * t57 - 0.2e1 * t56 * MDP(20) + MDP(25) * t106) * t57 + 0.2e1 * (-t77 * MDP(11) + pkin(1) * MDP(9)) * t91 + 0.2e1 * (-pkin(1) * MDP(10) - t77 * MDP(13) + t91 * MDP(5)) * t89; t89 * MDP(6) + t91 * MDP(7) + (-t89 * pkin(2) + t91 * qJ(3)) * MDP(12) + t61 * MDP(15) + t63 * MDP(16) + (-t76 * t67 + t74 * t68) * MDP(17) + (t61 * t74 + t63 * t76) * MDP(18) + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t91 + (-MDP(9) + t95) * t89) * pkin(6) - t93; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + (t74 ^ 2 + t76 ^ 2) * MDP(18) + MDP(23) + 0.2e1 * t74 * MDP(15) + 0.2e1 * t76 * MDP(16) + 0.2e1 * t101 + 0.2e1 * t100; (-t86 * t67 - t87 * t68) * MDP(17) + (-t61 * t87 + t63 * t86) * MDP(18) + (pkin(6) * MDP(14) + MDP(12)) * t89; -t87 * MDP(15) + t86 * MDP(16) + (-t74 * t87 + t76 * t86) * MDP(18) - t94 + t95; MDP(14) + (t86 ^ 2 + t87 ^ 2) * MDP(18); t57 * MDP(25) + t102 + t97 + t98 + t99; 0; 0; MDP(18); t93; -MDP(23) - t100 - t101; t94; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
