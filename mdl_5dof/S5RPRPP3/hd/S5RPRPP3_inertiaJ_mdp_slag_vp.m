% Calculate joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPP3_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:12:55
% DurationCPUTime: 0.15s
% Computational Cost: add. (235->77), mult. (399->97), div. (0->0), fcn. (376->4), ass. (0->27)
t65 = cos(qJ(3));
t64 = pkin(1) * MDP(7);
t50 = pkin(3) + qJ(5);
t63 = pkin(6) + qJ(2);
t48 = sin(pkin(7));
t49 = cos(pkin(7));
t52 = sin(qJ(3));
t40 = t52 * t48 - t65 * t49;
t61 = qJ(4) * t40;
t60 = t48 * MDP(5);
t59 = t49 * MDP(4);
t58 = MDP(17) + MDP(20);
t57 = (MDP(18) + MDP(22));
t45 = -t49 * pkin(2) - pkin(1);
t42 = t63 * t48;
t43 = t63 * t49;
t36 = t65 * t42 + t52 * t43;
t56 = -pkin(3) * MDP(18) + MDP(16);
t41 = t65 * t48 + t52 * t49;
t55 = -t41 * qJ(4) + t45;
t37 = -t52 * t42 + t65 * t43;
t53 = qJ(4) ^ 2;
t35 = t40 * pkin(3) + t55;
t34 = -t40 * pkin(4) + t37;
t33 = t41 * pkin(4) + t36;
t32 = t50 * t40 + t55;
t1 = [MDP(1) + (t35 ^ 2 + t36 ^ 2 + t37 ^ 2) * MDP(18) + (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) * MDP(22) + (0.2e1 * t59 - 0.2e1 * t60 + t64) * pkin(1) + 0.2e1 * (t45 * MDP(13) - t35 * MDP(16) + t32 * MDP(21)) * t40 + (0.2e1 * t45 * MDP(14) - 0.2e1 * t35 * MDP(17) - 0.2e1 * t32 * MDP(20) + MDP(8) * t41 - 0.2e1 * t40 * MDP(9)) * t41 + 0.2e1 * (t36 * t41 - t37 * t40) * MDP(15) + 0.2e1 * (t33 * t41 - t34 * t40) * MDP(19) + (MDP(7) * qJ(2) + 2 * MDP(6)) * (t48 ^ 2 + t49 ^ 2) * qJ(2); t35 * MDP(18) + t32 * MDP(22) - t59 + t60 - t64 + (MDP(14) - t58) * t41 + (MDP(13) - MDP(16) + MDP(21)) * t40; MDP(7) + t57; t41 * MDP(10) - t40 * MDP(11) + (-pkin(3) * t41 - t61) * MDP(15) + (-t50 * t41 - t61) * MDP(19) + t34 * MDP(20) - t33 * MDP(21) + (t34 * qJ(4) - t33 * t50) * MDP(22) + (qJ(4) * MDP(18) - MDP(14) + MDP(17)) * t37 + (-MDP(13) + t56) * t36; 0; MDP(12) - 0.2e1 * pkin(3) * MDP(16) + (pkin(3) ^ 2 + t53) * MDP(18) + 0.2e1 * t50 * MDP(21) + (t50 ^ 2 + t53) * MDP(22) + 0.2e1 * t58 * qJ(4); t36 * MDP(18) + t33 * MDP(22) + (MDP(15) + MDP(19)) * t41; 0; -t50 * MDP(22) - MDP(21) + t56; t57; -t40 * MDP(19) + t34 * MDP(22); 0; qJ(4) * MDP(22) + MDP(20); 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
