% Calculate joint inertia matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP3_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (118->50), mult. (207->70), div. (0->0), fcn. (116->4), ass. (0->27)
t51 = sin(qJ(2));
t40 = pkin(1) * t51 + pkin(6);
t50 = sin(qJ(3));
t48 = t50 ^ 2;
t52 = cos(qJ(3));
t64 = t52 ^ 2 + t48;
t66 = t64 * t40;
t70 = 2 * MDP(15);
t53 = cos(qJ(2));
t69 = t53 * pkin(1);
t41 = -pkin(2) - t69;
t68 = pkin(2) - t41;
t35 = -pkin(3) * t52 - qJ(4) * t50 - pkin(2);
t33 = t35 - t69;
t67 = -t33 - t35;
t65 = t64 * pkin(6);
t63 = MDP(17) * t50;
t62 = t50 * MDP(9) + t52 * MDP(10) + (-pkin(3) * t50 + qJ(4) * t52) * MDP(15);
t61 = 0.2e1 * t50 * t52 * MDP(8) + t48 * MDP(7) + MDP(4);
t60 = t64 * MDP(17);
t59 = -MDP(17) * pkin(3) - MDP(14);
t58 = MDP(12) * t52 - MDP(13) * t50;
t57 = -0.2e1 * MDP(14) * t52 - 0.2e1 * MDP(16) * t50;
t56 = (MDP(5) * t53 - MDP(6) * t51) * pkin(1);
t55 = (MDP(17) * qJ(4) - MDP(13) + MDP(16)) * t52 + (-MDP(12) + t59) * t50;
t42 = t50 * MDP(15);
t1 = [MDP(1) + t66 * t70 + t40 ^ 2 * t60 + (MDP(17) * t33 + t57) * t33 + t61 - 0.2e1 * t58 * t41 + 0.2e1 * t56; (t65 + t66) * MDP(15) + (pkin(6) * t66 + t33 * t35) * MDP(17) + t56 + (MDP(12) * t68 + MDP(14) * t67) * t52 + (-MDP(13) * t68 + MDP(16) * t67) * t50 + t61; t65 * t70 + pkin(6) ^ 2 * t60 + (MDP(17) * t35 + t57) * t35 + 0.2e1 * t58 * pkin(2) + t61; t40 * t55 + t62; pkin(6) * t55 + t62; MDP(11) + 0.2e1 * pkin(3) * MDP(14) + 0.2e1 * qJ(4) * MDP(16) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(17); t40 * t63 + t42; pkin(6) * t63 + t42; t59; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
