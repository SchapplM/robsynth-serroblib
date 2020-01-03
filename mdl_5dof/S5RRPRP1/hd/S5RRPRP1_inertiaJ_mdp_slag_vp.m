% Calculate joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:12
% EndTime: 2020-01-03 11:59:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (175->63), mult. (288->90), div. (0->0), fcn. (218->6), ass. (0->36)
t80 = 2 * MDP(15);
t65 = sin(qJ(2));
t79 = pkin(1) * t65;
t66 = cos(qJ(4));
t78 = t66 * pkin(4);
t67 = cos(qJ(2));
t56 = t67 * pkin(1) + pkin(2);
t62 = sin(pkin(8));
t63 = cos(pkin(8));
t45 = t63 * t56 - t62 * t79;
t42 = -pkin(3) - t45;
t55 = -t63 * pkin(2) - pkin(3);
t77 = t42 + t55;
t46 = t62 * t56 + t63 * t79;
t64 = sin(qJ(4));
t76 = t64 * MDP(10) + t66 * MDP(11);
t75 = t64 * MDP(14);
t74 = t64 * MDP(15);
t73 = t66 * MDP(13);
t61 = t64 ^ 2;
t72 = 0.2e1 * t64 * t66 * MDP(9) + t61 * MDP(8) + MDP(4);
t71 = -MDP(13) * t64 - MDP(14) * t66;
t70 = (t67 * MDP(5) - t65 * MDP(6)) * pkin(1);
t69 = -0.2e1 * t73 + 0.2e1 * t75;
t60 = t66 * qJ(5);
t54 = t62 * pkin(2) + pkin(7);
t49 = t55 - t78;
t48 = t66 * t54 + t60;
t47 = (-qJ(5) - t54) * t64;
t44 = t48 * t66;
t43 = pkin(7) + t46;
t41 = t42 - t78;
t40 = t66 * t43 + t60;
t39 = (-qJ(5) - t43) * t64;
t38 = t40 * t66;
t1 = [MDP(1) + (t45 ^ 2 + t46 ^ 2) * MDP(7) + (-t39 * t64 + t38) * t80 + (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) * MDP(16) + t42 * t69 + 0.2e1 * t70 + t72; (t38 + t44) * MDP(15) + (t39 * t47 + t40 * t48 + t41 * t49) * MDP(16) - t77 * t73 + (t45 * t63 + t46 * t62) * MDP(7) * pkin(2) + t70 + (t77 * MDP(14) + (-t39 - t47) * MDP(15)) * t64 + t72; (-t47 * t64 + t44) * t80 + (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) * MDP(16) + (t62 ^ 2 + t63 ^ 2) * MDP(7) * pkin(2) ^ 2 + t55 * t69 + t72; (t39 * t66 + t40 * t64) * MDP(16); (t47 * t66 + t48 * t64) * MDP(16); MDP(7) + (t66 ^ 2 + t61) * MDP(16); t71 * t43 + (t39 * MDP(16) - t74) * pkin(4) + t76; t71 * t54 + (MDP(16) * t47 - t74) * pkin(4) + t76; -t75 + (MDP(16) * pkin(4) + MDP(13)) * t66; pkin(4) ^ 2 * MDP(16) + MDP(12); t41 * MDP(16); t49 * MDP(16); 0; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
