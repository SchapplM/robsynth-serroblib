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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR5_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:58
% EndTime: 2019-12-31 19:29:59
% DurationCPUTime: 0.18s
% Computational Cost: add. (258->77), mult. (452->109), div. (0->0), fcn. (448->6), ass. (0->35)
t71 = sin(pkin(8));
t72 = cos(pkin(8));
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t59 = t71 * t74 - t72 * t76;
t60 = t71 * t76 + t72 * t74;
t70 = -pkin(2) * t76 - pkin(1);
t79 = qJ(4) * t60 - t70;
t92 = 0.2e1 * (-pkin(3) - pkin(4)) * t59 + 0.2e1 * t79;
t91 = 0.2e1 * t76;
t90 = -qJ(3) - pkin(6);
t64 = t90 * t76;
t82 = t90 * t74;
t55 = -t72 * t64 + t71 * t82;
t49 = pkin(3) * t59 - t79;
t89 = MDP(16) * t49;
t68 = pkin(2) * t72 + pkin(3);
t65 = -pkin(4) - t68;
t66 = pkin(2) * t71 + qJ(4);
t73 = sin(qJ(5));
t75 = cos(qJ(5));
t88 = MDP(22) * (-t65 * t75 + t66 * t73);
t87 = MDP(23) * (t65 * t73 + t66 * t75);
t50 = -t75 * t59 + t60 * t73;
t86 = t50 * MDP(22);
t85 = t59 * MDP(13);
t84 = t60 * MDP(15);
t53 = -t64 * t71 - t72 * t82;
t83 = t53 ^ 2 + t55 ^ 2;
t80 = MDP(22) * t75 - MDP(23) * t73;
t46 = -pkin(7) * t60 + t53;
t47 = pkin(7) * t59 + t55;
t51 = t59 * t73 + t60 * t75;
t78 = t51 * MDP(19) - t50 * MDP(20) - (-t46 * t75 + t47 * t73) * MDP(22) - (t46 * t73 + t47 * t75) * MDP(23);
t1 = [MDP(1) + pkin(1) * MDP(9) * t91 + (t70 ^ 2 + t83) * MDP(12) + t83 * MDP(16) + t86 * t92 + (-0.2e1 * t84 + 0.2e1 * t85 + t89) * t49 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t74 + MDP(5) * t91) * t74 + (MDP(17) * t51 - 0.2e1 * t50 * MDP(18) + MDP(23) * t92) * t51 + 0.2e1 * (MDP(11) + MDP(14)) * (t53 * t60 - t55 * t59); t74 * MDP(6) + t76 * MDP(7) - t53 * MDP(13) + (-t59 * t66 - t60 * t68) * MDP(14) + t55 * MDP(15) + (-t53 * t68 + t55 * t66) * MDP(16) + (-MDP(10) * t76 - MDP(9) * t74) * pkin(6) + ((-t59 * t71 - t60 * t72) * MDP(11) + (-t53 * t72 + t55 * t71) * MDP(12)) * pkin(2) - t78; MDP(8) + (t66 ^ 2 + t68 ^ 2) * MDP(16) + MDP(21) + (t71 ^ 2 + t72 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * MDP(13) * t68 + 0.2e1 * MDP(15) * t66 + 0.2e1 * t88 + 0.2e1 * t87; MDP(12) * t70 - t51 * MDP(23) - t84 + t85 - t86 + t89; 0; MDP(12) + MDP(16); t60 * MDP(14) + t53 * MDP(16); -MDP(16) * t68 - MDP(13) - t80; 0; MDP(16); t78; -MDP(21) - t87 - t88; 0; t80; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
