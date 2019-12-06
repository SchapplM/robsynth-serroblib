% Calculate joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:01
% EndTime: 2019-12-05 18:46:01
% DurationCPUTime: 0.24s
% Computational Cost: add. (415->88), mult. (747->132), div. (0->0), fcn. (790->6), ass. (0->45)
t78 = cos(qJ(2));
t71 = -t78 * pkin(2) - pkin(1);
t97 = 0.2e1 * t71;
t96 = 0.2e1 * t78;
t95 = pkin(6) + pkin(7);
t74 = sin(qJ(3));
t94 = pkin(2) * t74;
t73 = sin(qJ(4));
t93 = t73 * pkin(3);
t92 = MDP(26) * pkin(3);
t76 = cos(qJ(4));
t72 = t76 * pkin(3);
t69 = t72 + pkin(4);
t91 = MDP(26) * t69;
t77 = cos(qJ(3));
t70 = t77 * pkin(2) + pkin(3);
t68 = t76 * t70;
t58 = -t73 * t94 + t68;
t90 = t58 * MDP(23);
t59 = t73 * t70 + t76 * t94;
t89 = t59 * MDP(24);
t88 = t76 * MDP(23);
t87 = t77 * MDP(16);
t86 = MDP(15) + MDP(22);
t75 = sin(qJ(2));
t64 = t74 * t78 + t77 * t75;
t66 = t95 * t75;
t67 = t95 * t78;
t84 = -t77 * t66 - t74 * t67;
t48 = -t64 * pkin(8) + t84;
t63 = t74 * t75 - t77 * t78;
t81 = t74 * t66 - t77 * t67;
t49 = -t63 * pkin(8) - t81;
t85 = t76 * t48 - t73 * t49;
t52 = t76 * t63 + t73 * t64;
t53 = -t73 * t63 + t76 * t64;
t82 = -t73 * t48 - t76 * t49;
t83 = t53 * MDP(20) - t52 * MDP(21) + t85 * MDP(23) + t82 * MDP(24);
t80 = t63 * pkin(3) + t71;
t79 = t64 * MDP(13) - t63 * MDP(14) + t84 * MDP(16) + t81 * MDP(17) + t83;
t57 = pkin(4) + t58;
t46 = t52 * pkin(4) + t80;
t43 = -t52 * qJ(5) - t82;
t42 = -t53 * qJ(5) + t85;
t1 = [MDP(1) + pkin(1) * MDP(9) * t96 + t63 * MDP(16) * t97 + t53 ^ 2 * MDP(18) - 0.2e1 * t53 * t52 * MDP(19) + 0.2e1 * (-t42 * t53 - t43 * t52) * MDP(25) + (t42 ^ 2 + t43 ^ 2 + t46 ^ 2) * MDP(26) + 0.2e1 * (t52 * MDP(23) + t53 * MDP(24)) * t80 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t75 + MDP(5) * t96) * t75 + (MDP(11) * t64 - 0.2e1 * t63 * MDP(12) + MDP(17) * t97) * t64; t75 * MDP(6) + t78 * MDP(7) + (-t59 * t52 - t57 * t53) * MDP(25) + (t42 * t57 + t43 * t59) * MDP(26) + (-t78 * MDP(10) - t75 * MDP(9)) * pkin(6) + t79; MDP(8) + (t57 ^ 2 + t59 ^ 2) * MDP(26) + 0.2e1 * (-t74 * MDP(17) + t87) * pkin(2) + 0.2e1 * t90 - 0.2e1 * t89 + t86; (-t52 * t93 - t69 * t53) * MDP(25) + (t42 * t69 + t43 * t93) * MDP(26) + t79; (t68 + t72) * MDP(23) + t57 * t91 + ((-pkin(3) - t70) * MDP(24) + t59 * t92) * t73 + (t87 + (-MDP(23) * t73 - MDP(24) * t76 - MDP(17)) * t74) * pkin(2) + t86; t69 ^ 2 * MDP(26) + (0.2e1 * t88 + (t73 * t92 - 0.2e1 * MDP(24)) * t73) * pkin(3) + t86; (-t53 * MDP(25) + t42 * MDP(26)) * pkin(4) + t83; t57 * pkin(4) * MDP(26) + MDP(22) - t89 + t90; pkin(4) * t91 + MDP(22) + (-t73 * MDP(24) + t88) * pkin(3); MDP(26) * pkin(4) ^ 2 + MDP(22); t46 * MDP(26); 0; 0; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
