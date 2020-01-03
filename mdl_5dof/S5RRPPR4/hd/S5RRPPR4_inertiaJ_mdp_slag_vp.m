% Calculate joint inertia matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRPPR4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (186->70), mult. (235->82), div. (0->0), fcn. (158->6), ass. (0->36)
t68 = cos(qJ(5));
t79 = t68 * MDP(18);
t66 = sin(qJ(5));
t80 = t66 * MDP(19);
t73 = t79 - t80;
t87 = 2 * MDP(7);
t86 = 2 * qJ(3);
t85 = 2 * MDP(10);
t84 = 2 * MDP(11);
t69 = cos(qJ(2));
t58 = t69 * pkin(1) + pkin(2);
t54 = -pkin(3) - t58;
t67 = sin(qJ(2));
t63 = t67 * pkin(1);
t55 = t63 + qJ(3);
t64 = sin(pkin(8));
t65 = cos(pkin(8));
t43 = t64 * t54 + t65 * t55;
t70 = -pkin(2) - pkin(3);
t48 = t65 * qJ(3) + t64 * t70;
t82 = t67 * MDP(6);
t81 = MDP(12) * t64;
t78 = MDP(4) + (MDP(13) * t66 + 0.2e1 * MDP(14) * t68) * t66;
t77 = t64 * MDP(11) + t65 * t80 - MDP(7);
t50 = t65 * t54;
t41 = t64 * t55 - t50;
t57 = t65 * t70;
t46 = t64 * qJ(3) - t57;
t76 = -MDP(10) - t79;
t75 = pkin(2) * t87 + t78;
t74 = -t66 * MDP(15) - t68 * MDP(16);
t72 = -MDP(18) * t66 - MDP(19) * t68;
t71 = 0.2e1 * t73;
t44 = pkin(4) + t46;
t39 = pkin(4) + t41;
t1 = [MDP(1) + (t55 ^ 2 + t58 ^ 2) * MDP(9) + (t41 ^ 2 + t43 ^ 2) * MDP(12) + t39 * t71 + 0.2e1 * (t69 * MDP(5) - t82) * pkin(1) + t41 * t85 + t43 * t84 + t58 * t87 + 0.2e1 * t55 * MDP(8) + t78; (t86 + t63) * MDP(8) + (t58 * pkin(2) + t55 * qJ(3)) * MDP(9) + (-t50 - t57 + (qJ(3) + t55) * t64) * MDP(10) + (t48 + t43) * MDP(11) + (t41 * t46 + t43 * t48) * MDP(12) + (-t82 + (MDP(5) + MDP(7)) * t69) * pkin(1) + t75 + t73 * (t39 + t44); MDP(8) * t86 + (pkin(2) ^ 2 + (qJ(3) ^ 2)) * MDP(9) + (t46 ^ 2 + t48 ^ 2) * MDP(12) + t44 * t71 + t46 * t85 + t48 * t84 + t75; t43 * t81 - t58 * MDP(9) + (-t41 * MDP(12) + t76) * t65 + t77; t48 * t81 - pkin(2) * MDP(9) + (-t46 * MDP(12) + t76) * t65 + t77; MDP(9) + (t64 ^ 2 + t65 ^ 2) * MDP(12); 0; 0; 0; MDP(12); t72 * (-pkin(7) + t43) + t74; t72 * (-pkin(7) + t48) + t74; t72 * t64; t73; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
