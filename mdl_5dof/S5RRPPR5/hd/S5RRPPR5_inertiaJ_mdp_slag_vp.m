% Calculate joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR5_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:32
% EndTime: 2021-01-15 19:36:33
% DurationCPUTime: 0.23s
% Computational Cost: add. (278->82), mult. (490->114), div. (0->0), fcn. (480->6), ass. (0->38)
t95 = 2 * MDP(11);
t94 = MDP(14) + MDP(18);
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t73 = sin(qJ(2));
t75 = cos(qJ(2));
t60 = t70 * t73 - t71 * t75;
t61 = t70 * t75 + t71 * t73;
t69 = -t75 * pkin(2) - pkin(1);
t78 = t61 * qJ(4) - t69;
t93 = 0.2e1 * (-pkin(3) - pkin(4)) * t60 + 0.2e1 * t78;
t92 = 0.2e1 * t75;
t91 = 2 * MDP(15);
t90 = -qJ(3) - pkin(6);
t50 = t60 * pkin(3) - t78;
t89 = t50 * MDP(18);
t72 = sin(qJ(5));
t74 = cos(qJ(5));
t51 = -t74 * t60 + t72 * t61;
t88 = t51 * MDP(24);
t67 = t71 * pkin(2) + pkin(3);
t64 = -pkin(4) - t67;
t65 = t70 * pkin(2) + qJ(4);
t87 = (-t74 * t64 + t72 * t65) * MDP(24);
t86 = (t72 * t64 + t74 * t65) * MDP(25);
t85 = t69 * MDP(14);
t84 = MDP(12) - MDP(17);
t82 = t90 * t73;
t63 = t90 * t75;
t54 = -t70 * t63 - t71 * t82;
t80 = -t67 * MDP(18) - MDP(15);
t79 = t74 * MDP(24) - t72 * MDP(25);
t56 = -t71 * t63 + t70 * t82;
t47 = -t61 * pkin(7) + t54;
t48 = t60 * pkin(7) + t56;
t52 = t72 * t60 + t74 * t61;
t77 = t52 * MDP(21) - t51 * MDP(22) - (-t74 * t47 + t72 * t48) * MDP(24) - (t72 * t47 + t74 * t48) * MDP(25);
t1 = [MDP(1) + pkin(1) * MDP(9) * t92 + t88 * t93 + (0.2e1 * t61 * MDP(12) + t60 * t95 + t85) * t69 + (-0.2e1 * t61 * MDP(17) + t60 * t91 + t89) * t50 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t73 + MDP(5) * t92) * t73 + (MDP(19) * t52 - 0.2e1 * t51 * MDP(20) + MDP(25) * t93) * t52 + t94 * (t54 ^ 2 + t56 ^ 2) + 0.2e1 * (MDP(13) + MDP(16)) * (t54 * t61 - t56 * t60); t73 * MDP(6) + t75 * MDP(7) + (-t65 * t60 - t67 * t61) * MDP(16) + (MDP(18) * t65 - t84) * t56 + (-MDP(11) + t80) * t54 + (-t75 * MDP(10) - t73 * MDP(9)) * pkin(6) + ((-t60 * t70 - t61 * t71) * MDP(13) + (-t54 * t71 + t56 * t70) * MDP(14)) * pkin(2) - t77; MDP(8) + (t65 ^ 2 + t67 ^ 2) * MDP(18) + MDP(23) + t67 * t91 + 0.2e1 * t65 * MDP(17) + 0.2e1 * t87 + 0.2e1 * t86 + (t71 * t95 - 0.2e1 * t70 * MDP(12) + (t70 ^ 2 + t71 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t85 + t89 - t88 - t52 * MDP(25) + t84 * t61 + (MDP(11) + MDP(15)) * t60; 0; t94; t61 * MDP(16) + t54 * MDP(18); -t79 + t80; 0; MDP(18); t77; -MDP(23) - t86 - t87; 0; t79; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
