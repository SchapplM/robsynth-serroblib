% Calculate joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR8_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:24
% DurationCPUTime: 0.12s
% Computational Cost: add. (86->50), mult. (154->74), div. (0->0), fcn. (107->4), ass. (0->21)
t43 = cos(qJ(2));
t41 = sin(qJ(2));
t49 = t41 * qJ(3) + pkin(1);
t54 = pkin(2) + pkin(3);
t55 = 0.2e1 * t54 * t43 + 0.2e1 * t49;
t53 = pkin(5) - pkin(6);
t38 = t41 ^ 2;
t52 = t43 ^ 2 + t38;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t51 = (t40 * qJ(3) + t42 * t54) * MDP(20);
t50 = (t42 * qJ(3) - t40 * t54) * MDP(21);
t48 = -pkin(2) * MDP(14) - MDP(11);
t47 = t42 * MDP(20) - t40 * MDP(21);
t31 = t41 * t40 + t43 * t42;
t32 = -t43 * t40 + t41 * t42;
t36 = t53 * t41;
t37 = t53 * t43;
t46 = t32 * MDP(17) - t31 * MDP(18) - (-t42 * t36 + t40 * t37) * MDP(20) - (t40 * t36 + t42 * t37) * MDP(21);
t35 = -t43 * pkin(2) - t49;
t1 = [MDP(1) + t38 * MDP(4) + (t52 * pkin(5) ^ 2 + t35 ^ 2) * MDP(14) + t31 * MDP(20) * t55 + 0.2e1 * t52 * MDP(12) * pkin(5) + (MDP(15) * t32 - 0.2e1 * t31 * MDP(16) + MDP(21) * t55) * t32 + 0.2e1 * (-t35 * MDP(11) + pkin(1) * MDP(9)) * t43 + 0.2e1 * (-pkin(1) * MDP(10) - t35 * MDP(13) + t43 * MDP(5)) * t41; t41 * MDP(6) + t43 * MDP(7) + (-t41 * pkin(2) + t43 * qJ(3)) * MDP(12) + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t43 + (-MDP(9) + t48) * t41) * pkin(5) - t46; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + MDP(19) + 0.2e1 * t51 + 0.2e1 * t50; (MDP(14) * pkin(5) + MDP(12)) * t41; -t47 + t48; MDP(14); t46; -MDP(19) - t50 - t51; t47; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
