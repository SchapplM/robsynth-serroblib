% Calculate joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR11_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:54
% EndTime: 2019-12-31 18:05:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (95->60), mult. (166->77), div. (0->0), fcn. (110->4), ass. (0->25)
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t60 = MDP(22) * t41 + MDP(23) * t43;
t44 = cos(qJ(4));
t59 = t44 * MDP(16) + MDP(8);
t38 = -pkin(6) + qJ(2);
t57 = t38 * t41;
t56 = t38 * t43;
t55 = t41 * t43;
t39 = pkin(1) + qJ(3);
t51 = MDP(18) * t55;
t50 = MDP(19) * t43 - MDP(20) * t41;
t49 = t43 * MDP(22) - t41 * MDP(23);
t47 = MDP(15) + t49;
t46 = t41 * MDP(19) + t43 * MDP(20) - pkin(7) * t60;
t45 = (qJ(2) ^ 2);
t42 = sin(qJ(4));
t37 = t44 ^ 2;
t36 = t43 ^ 2;
t35 = t42 ^ 2;
t34 = t41 ^ 2;
t33 = pkin(4) * t42 - pkin(7) * t44 + t39;
t32 = t33 * t41 + t42 * t56;
t31 = t33 * t43 - t42 * t57;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t45) * MDP(6)) + (t39 ^ 2 + t45) * MDP(9) + t35 * MDP(21) + (t36 * MDP(17) + MDP(10) - 0.2e1 * t51) * t37 + 0.2e1 * (t31 * t42 - t37 * t57) * MDP(22) + 0.2e1 * (-t32 * t42 - t37 * t56) * MDP(23) + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (-MDP(11) + t50) * t42 * t44 + 0.2e1 * (t42 * MDP(15) + t59) * t39; -(pkin(1) * MDP(6)) - MDP(9) * t39 - t47 * t42 + MDP(4) - t59; MDP(6) + MDP(9); qJ(2) * MDP(9) + MDP(7) + t60 * (-t35 - t37); 0; MDP(9); (-t38 * MDP(16) - MDP(13) + t46) * t42 + (MDP(12) + t38 * MDP(15) + MDP(17) * t55 + (-t34 + t36) * MDP(18) + (-pkin(4) * t41 + t56) * MDP(22) + (-pkin(4) * t43 - t57) * MDP(23)) * t44; 0; -MDP(16) * t42 + t47 * t44; MDP(17) * t34 + 0.2e1 * pkin(4) * t49 + MDP(14) + 0.2e1 * t51; MDP(21) * t42 + t31 * MDP(22) - t32 * MDP(23) + t50 * t44; -t49; -t60 * t42; t46; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
