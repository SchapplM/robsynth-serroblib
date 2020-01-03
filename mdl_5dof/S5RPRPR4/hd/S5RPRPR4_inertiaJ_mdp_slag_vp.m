% Calculate joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR4_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:18
% EndTime: 2020-01-03 11:39:19
% DurationCPUTime: 0.18s
% Computational Cost: add. (232->60), mult. (412->99), div. (0->0), fcn. (431->8), ass. (0->32)
t65 = sin(pkin(9));
t67 = cos(pkin(9));
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t59 = -t65 * t70 + t67 * t72;
t68 = cos(pkin(8));
t64 = -t68 * pkin(1) - pkin(2);
t61 = -t72 * pkin(3) + t64;
t82 = -0.2e1 * t59 * pkin(4) + 0.2e1 * t61;
t81 = pkin(3) * t65;
t60 = t65 * t72 + t67 * t70;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t48 = -t71 * t59 + t69 * t60;
t44 = t48 * MDP(19);
t49 = t69 * t59 + t71 * t60;
t80 = -t49 * MDP(20) - t44;
t66 = sin(pkin(8));
t62 = t66 * pkin(1) + pkin(6);
t79 = qJ(4) + t62;
t57 = t79 * t70;
t58 = t79 * t72;
t43 = -t65 * t57 + t67 * t58;
t63 = t67 * pkin(3) + pkin(4);
t78 = (t71 * t63 - t69 * t81) * MDP(19);
t77 = (-t69 * t63 - t71 * t81) * MDP(20);
t76 = t72 * MDP(10);
t42 = -t67 * t57 - t65 * t58;
t40 = -t60 * pkin(7) + t42;
t41 = t59 * pkin(7) + t43;
t75 = t49 * MDP(16) - t48 * MDP(17) + (t71 * t40 - t69 * t41) * MDP(19) + (-t69 * t40 - t71 * t41) * MDP(20);
t1 = [MDP(1) - 0.2e1 * t64 * t76 + 0.2e1 * (-t42 * t60 + t43 * t59) * MDP(12) + (t42 ^ 2 + t43 ^ 2 + t61 ^ 2) * MDP(13) + t44 * t82 + (t66 ^ 2 + t68 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * t64 * MDP(11) + MDP(5) * t70 + 0.2e1 * t72 * MDP(6)) * t70 + (MDP(14) * t49 - 0.2e1 * t48 * MDP(15) + MDP(20) * t82) * t49; (t42 * t59 + t43 * t60) * MDP(13); MDP(4) + (t59 ^ 2 + t60 ^ 2) * MDP(13); t70 * MDP(7) + t72 * MDP(8) + (-t70 * MDP(10) - t72 * MDP(11)) * t62 + ((t59 * t65 - t60 * t67) * MDP(12) + (t42 * t67 + t43 * t65) * MDP(13)) * pkin(3) + t75; t76 - t70 * MDP(11) + (t59 * t67 + t60 * t65) * MDP(13) * pkin(3) + t80; MDP(9) + MDP(18) + (t65 ^ 2 + t67 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t78 + 0.2e1 * t77; t61 * MDP(13) - t80; 0; 0; MDP(13); t75; t80; MDP(18) + t77 + t78; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
