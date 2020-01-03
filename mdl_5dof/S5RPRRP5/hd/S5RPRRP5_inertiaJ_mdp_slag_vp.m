% Calculate joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP5_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:58
% EndTime: 2019-12-31 18:40:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (186->59), mult. (314->82), div. (0->0), fcn. (214->6), ass. (0->35)
t62 = cos(pkin(8));
t52 = t62 * pkin(1) + pkin(2);
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t61 = sin(pkin(8));
t86 = pkin(1) * t61;
t47 = -t64 * t52 - t66 * t86;
t45 = pkin(7) - t47;
t63 = sin(qJ(4));
t59 = t63 ^ 2;
t65 = cos(qJ(4));
t81 = t65 ^ 2 + t59;
t83 = t81 * t45;
t87 = 2 * MDP(16);
t46 = t66 * t52 - t64 * t86;
t44 = -pkin(3) - t46;
t85 = pkin(3) - t44;
t49 = -t65 * pkin(4) - t63 * qJ(5) - pkin(3);
t40 = -t46 + t49;
t84 = -t40 - t49;
t82 = t81 * pkin(7);
t80 = t46 * MDP(6);
t79 = t47 * MDP(7);
t78 = MDP(18) * t63;
t77 = t63 * MDP(10) + t65 * MDP(11) + (-t63 * pkin(4) + t65 * qJ(5)) * MDP(16);
t76 = 0.2e1 * t63 * t65 * MDP(9) + t59 * MDP(8) + MDP(5);
t75 = t81 * MDP(18);
t74 = -pkin(4) * MDP(18) - MDP(15);
t73 = MDP(13) - t74;
t72 = MDP(18) * qJ(5) - MDP(14) + MDP(17);
t71 = t65 * MDP(13) - t63 * MDP(14);
t70 = -0.2e1 * t65 * MDP(15) - 0.2e1 * t63 * MDP(17);
t69 = -t63 * t73 + t65 * t72;
t53 = t63 * MDP(16);
t1 = [MDP(1) + (t61 ^ 2 + t62 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t71 * t44 + t45 ^ 2 * t75 + (t40 * MDP(18) + t70) * t40 + 0.2e1 * t80 + 0.2e1 * t79 + t83 * t87 + t76; 0; MDP(4) + t75; t80 + t79 + (t82 + t83) * MDP(16) + (pkin(7) * t83 + t40 * t49) * MDP(18) + (MDP(13) * t85 + MDP(15) * t84) * t65 + (-MDP(14) * t85 + MDP(17) * t84) * t63 + t76; 0; t82 * t87 + pkin(7) ^ 2 * t75 + (MDP(18) * t49 + t70) * t49 + 0.2e1 * t71 * pkin(3) + t76; t45 * t69 + t77; t63 * t72 + t65 * t73; pkin(7) * t69 + t77; MDP(12) + 0.2e1 * pkin(4) * MDP(15) + 0.2e1 * qJ(5) * MDP(17) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(18); t45 * t78 + t53; -t65 * MDP(18); pkin(7) * t78 + t53; t74; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
