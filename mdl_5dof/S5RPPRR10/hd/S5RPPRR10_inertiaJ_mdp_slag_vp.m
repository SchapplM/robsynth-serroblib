% Calculate joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR10_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:18
% EndTime: 2019-12-31 18:04:19
% DurationCPUTime: 0.21s
% Computational Cost: add. (172->60), mult. (315->85), div. (0->0), fcn. (314->6), ass. (0->31)
t91 = (MDP(7) + MDP(11));
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t59 = -t71 * pkin(2) - t70 * qJ(3) - pkin(1);
t53 = t71 * pkin(3) - t59;
t73 = sin(qJ(4));
t75 = cos(qJ(4));
t57 = t70 * t73 + t71 * t75;
t90 = 0.2e1 * t57 * pkin(4) + 0.2e1 * t53;
t89 = 0.2e1 * t53;
t88 = -0.2e1 * t70;
t87 = pkin(1) * MDP(7);
t86 = -pkin(6) + qJ(2);
t72 = sin(qJ(5));
t74 = cos(qJ(5));
t85 = (-t72 * t73 + t74 * t75) * MDP(24) + (-t72 * t75 - t74 * t73) * MDP(25);
t58 = t70 * t75 - t71 * t73;
t49 = t74 * t57 + t72 * t58;
t83 = t49 * MDP(24);
t82 = t57 * MDP(17);
t81 = t59 * MDP(11);
t61 = t86 * t70;
t62 = t86 * t71;
t80 = t75 * t61 - t73 * t62;
t45 = -t58 * pkin(7) + t80;
t78 = -t73 * t61 - t75 * t62;
t46 = -t57 * pkin(7) - t78;
t50 = -t72 * t57 + t74 * t58;
t79 = t50 * MDP(21) - t49 * MDP(22) + (t74 * t45 - t72 * t46) * MDP(24) + (-t72 * t45 - t74 * t46) * MDP(25);
t77 = (MDP(24) * t74 - MDP(25) * t72) * pkin(4);
t1 = [MDP(1) + t82 * t89 + t83 * t90 + (MDP(10) * t88 - 0.2e1 * t71 * MDP(8) + t81) * t59 + (0.2e1 * t71 * MDP(4) + MDP(5) * t88 + t87) * pkin(1) + (MDP(12) * t58 - 0.2e1 * t57 * MDP(13) + MDP(18) * t89) * t58 + (MDP(19) * t50 - 0.2e1 * t49 * MDP(20) + MDP(25) * t90) * t50 + (t91 * qJ(2) + 2 * MDP(6) + 2 * MDP(9)) * (t70 ^ 2 + t71 ^ 2) * qJ(2); t81 - t82 - t58 * MDP(18) - t83 - t50 * MDP(25) - t87 + (-MDP(4) - MDP(8)) * t71 + (-MDP(10) + MDP(5)) * t70; t91; (MDP(11) * qJ(2) + MDP(9)) * t70; 0; MDP(11); t58 * MDP(14) - t57 * MDP(15) + MDP(17) * t80 + MDP(18) * t78 + t79; 0; t75 * MDP(17) - t73 * MDP(18) + t85; MDP(16) + MDP(23) + 0.2e1 * t77; t79; 0; t85; MDP(23) + t77; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
