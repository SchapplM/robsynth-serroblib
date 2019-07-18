% Calculate joint inertia matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPR2_inertiaJ_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (84->45), mult. (104->51), div. (0->0), fcn. (55->4), ass. (0->23)
t47 = 2 * MDP(7);
t46 = 2 * qJ(3);
t33 = cos(qJ(2));
t26 = t33 * pkin(1) + pkin(2);
t24 = -pkin(3) - t26;
t34 = -pkin(2) - pkin(3);
t45 = t24 + t34;
t31 = sin(qJ(2));
t44 = t31 * MDP(6);
t29 = t31 * pkin(1);
t25 = t29 + qJ(3);
t43 = qJ(3) + t25;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t42 = (-t32 * t24 + t30 * t25) * MDP(11);
t41 = (t30 * t24 + t32 * t25) * MDP(12);
t40 = (t30 * qJ(3) - t32 * t34) * MDP(11);
t39 = (t32 * qJ(3) + t30 * t34) * MDP(12);
t38 = MDP(10) + MDP(4);
t37 = t32 * MDP(11) - t30 * MDP(12);
t36 = pkin(2) * t47 + t38;
t35 = -MDP(7) - t37;
t1 = [MDP(1) + (t25 ^ 2 + t26 ^ 2) * MDP(9) + 0.2e1 * (t33 * MDP(5) - t44) * pkin(1) + 0.2e1 * t42 + 0.2e1 * t41 + t26 * t47 + 0.2e1 * t25 * MDP(8) + t38; (t46 + t29) * MDP(8) + (t26 * pkin(2) + t25 * qJ(3)) * MDP(9) + (-t45 * MDP(11) + t43 * MDP(12)) * t32 + (t43 * MDP(11) + t45 * MDP(12)) * t30 + (-t44 + (MDP(5) + MDP(7)) * t33) * pkin(1) + t36; MDP(8) * t46 + (pkin(2) ^ 2 + (qJ(3) ^ 2)) * MDP(9) + 0.2e1 * t40 + 0.2e1 * t39 + t36; -t26 * MDP(9) + t35; -pkin(2) * MDP(9) + t35; MDP(9); -MDP(10) - t41 - t42; -MDP(10) - t39 - t40; t37; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq  = res;
