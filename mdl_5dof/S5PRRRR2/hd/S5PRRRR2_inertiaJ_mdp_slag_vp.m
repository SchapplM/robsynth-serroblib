% Calculate joint inertia matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:19
% EndTime: 2019-07-18 13:30:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (86->35), mult. (172->40), div. (0->0), fcn. (118->6), ass. (0->24)
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t44 = t37 * MDP(16) - t34 * MDP(17);
t38 = cos(qJ(4));
t42 = MDP(9) + t44;
t35 = sin(qJ(4));
t49 = t35 * MDP(10);
t56 = (t42 * t38 - t49) * pkin(3);
t36 = sin(qJ(3));
t55 = pkin(2) * t36;
t54 = t38 * pkin(3);
t53 = t34 * MDP(13) + t37 * MDP(14);
t39 = cos(qJ(3));
t52 = t39 * MDP(6);
t30 = t39 * pkin(2) + pkin(3);
t26 = -t35 * t30 - t38 * t55;
t51 = t26 * MDP(10);
t47 = MDP(8) + (MDP(11) * t34 + 0.2e1 * MDP(12) * t37) * t34;
t45 = MDP(5) + t47;
t43 = -MDP(16) * t34 - MDP(17) * t37;
t27 = t38 * t30;
t25 = t35 * t55 - t27;
t41 = t42 * t25;
t1 = [MDP(1); 0; 0.2e1 * t51 + MDP(2) + t45 + 0.2e1 * (-t36 * MDP(7) + t52) * pkin(2) - 0.2e1 * t41; 0; (t27 + t54) * MDP(9) + (-pkin(3) - t30) * t49 + (t52 + (-MDP(10) * t38 - MDP(9) * t35 - MDP(7)) * t36) * pkin(2) + t45 + t44 * (-t25 + t54); t45 + 0.2e1 * t56; 0; -t41 + t47 + t51; t56 + t47; t47; t44; t43 * (pkin(6) - t26) + t53; t43 * (t35 * pkin(3) + pkin(6)) + t53; t43 * pkin(6) + t53; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
