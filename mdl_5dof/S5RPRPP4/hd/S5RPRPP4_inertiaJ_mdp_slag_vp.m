% Calculate joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:51
% DurationCPUTime: 0.17s
% Computational Cost: add. (231->66), mult. (348->93), div. (0->0), fcn. (326->4), ass. (0->28)
t68 = sin(pkin(7));
t69 = cos(pkin(7));
t70 = cos(qJ(3));
t86 = sin(qJ(3));
t55 = -t68 * t86 + t69 * t70;
t56 = -t68 * t70 - t69 * t86;
t92 = t55 ^ 2 + t56 ^ 2;
t87 = -pkin(1) - pkin(6);
t74 = t86 * t87;
t59 = -t86 * qJ(4) + t74;
t75 = (-qJ(4) + t87) * t70;
t50 = t59 * t68 - t69 * t75;
t52 = t69 * t59 + t68 * t75;
t91 = t50 * t55 + t52 * t56;
t90 = MDP(14) + MDP(17);
t89 = MDP(15) + MDP(19);
t85 = (pkin(1) * MDP(6));
t66 = t86 * pkin(3) + qJ(2);
t82 = MDP(16) * t56;
t81 = MDP(18) * t55;
t48 = -pkin(4) * t56 - qJ(5) * t55 + t66;
t80 = MDP(19) * t48;
t79 = t50 ^ 2 + t52 ^ 2;
t60 = pkin(3) * t68 + qJ(5);
t64 = pkin(3) * t69 + pkin(4);
t73 = t55 * t64 - t56 * t60;
t72 = t55 * t69 - t56 * t68;
t1 = [MDP(1) + (t66 ^ 2 + t79) * MDP(15) + t79 * MDP(19) + (MDP(7) * t70 - 0.2e1 * t86 * MDP(8)) * t70 + (t80 - 0.2e1 * t81 - 0.2e1 * t82) * t48 + ((-2 * MDP(4) + t85) * pkin(1)) + (0.2e1 * t86 * MDP(12) + 0.2e1 * t70 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * t90 * t91; -t89 * t91 - t90 * t92 + MDP(4) - t85; t89 * t92 + MDP(6); -t86 * MDP(10) - MDP(13) * t74 - t50 * MDP(16) - t73 * MDP(17) + t52 * MDP(18) + (-t50 * t64 + t52 * t60) * MDP(19) + (t87 * MDP(12) + MDP(9)) * t70 + (-t72 * MDP(14) + (-t50 * t69 + t52 * t68) * MDP(15)) * pkin(3); t72 * MDP(15) * pkin(3) + t70 * MDP(12) - t86 * MDP(13) + t55 * MDP(16) - t56 * MDP(18) + t73 * MDP(19); MDP(11) + (t60 ^ 2 + t64 ^ 2) * MDP(19) + (t68 ^ 2 + t69 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * MDP(16) * t64 + 0.2e1 * MDP(18) * t60; MDP(15) * t66 + t80 - t81 - t82; 0; 0; t89; t55 * MDP(17) + MDP(19) * t50; -t55 * MDP(19); -MDP(19) * t64 - MDP(16); 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
