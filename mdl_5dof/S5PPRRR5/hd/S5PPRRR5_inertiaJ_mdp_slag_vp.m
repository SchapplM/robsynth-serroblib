% Calculate joint inertia matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR5_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:45
% DurationCPUTime: 0.09s
% Computational Cost: add. (52->25), mult. (101->32), div. (0->0), fcn. (85->6), ass. (0->15)
t31 = sin(qJ(5));
t34 = cos(qJ(5));
t40 = t34 * MDP(14) - t31 * MDP(15);
t44 = t31 * MDP(11) + t34 * MDP(12);
t41 = MDP(6) + (0.2e1 * MDP(10) * t34 + MDP(9) * t31) * t31;
t39 = -MDP(14) * t31 - MDP(15) * t34;
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t35 = cos(qJ(4));
t36 = cos(qJ(3));
t24 = t32 * t36 + t35 * t33;
t38 = -t24 * MDP(8) + (-MDP(7) - t40) * (t32 * t33 - t35 * t36);
t37 = (t35 * MDP(7) - t32 * MDP(8)) * pkin(3);
t27 = -t35 * pkin(3) - pkin(4);
t1 = [MDP(1) + MDP(2); 0; MDP(2); 0; t36 * MDP(4) - t33 * MDP(5) + t38; -0.2e1 * t27 * t40 + MDP(3) + 0.2e1 * t37 + t41; 0; t38; t37 + t41 + t40 * (pkin(4) - t27); 0.2e1 * pkin(4) * t40 + t41; -t40; t39 * t24; t39 * (t32 * pkin(3) + pkin(7)) + t44; t39 * pkin(7) + t44; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
