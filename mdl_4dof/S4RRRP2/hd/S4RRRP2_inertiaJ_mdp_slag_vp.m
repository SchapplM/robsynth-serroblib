% Calculate joint inertia matrix for
% S4RRRP2
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
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:30
% EndTime: 2021-01-15 11:04:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (124->57), mult. (204->71), div. (0->0), fcn. (133->4), ass. (0->31)
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t50 = -t44 * MDP(14) + t42 * MDP(15);
t62 = 2 * MDP(16);
t45 = cos(qJ(2));
t61 = t45 * pkin(1);
t36 = -pkin(2) - t61;
t60 = pkin(2) - t36;
t59 = -qJ(4) - pkin(6);
t37 = -t44 * pkin(3) - pkin(2);
t31 = t37 - t61;
t58 = t31 + t37;
t57 = t44 * MDP(10) + t42 * MDP(9);
t43 = sin(qJ(2));
t35 = t43 * pkin(1) + pkin(6);
t56 = qJ(4) + t35;
t55 = t31 * MDP(17);
t54 = t37 * MDP(17);
t53 = t42 * MDP(16);
t51 = MDP(4) + (MDP(7) * t42 + 0.2e1 * MDP(8) * t44) * t42;
t49 = t44 * MDP(12) - t42 * MDP(13);
t48 = -MDP(12) * t42 - MDP(13) * t44;
t47 = 0.2e1 * t50;
t46 = (t45 * MDP(5) - t43 * MDP(6)) * pkin(1);
t33 = t59 * t44;
t32 = t59 * t42;
t30 = t33 * t44;
t29 = t56 * t44;
t28 = t56 * t42;
t27 = t29 * t44;
t1 = [MDP(1) + (t28 * t42 + t27) * t62 + (t28 ^ 2 + t29 ^ 2) * MDP(17) + (t47 + t55) * t31 + t51 - 0.2e1 * t49 * t36 + 0.2e1 * t46; (t27 - t30) * MDP(16) + (-t28 * t32 - t29 * t33 + t31 * t37) * MDP(17) + t46 + (t60 * MDP(12) - t58 * MDP(14)) * t44 + (-t60 * MDP(13) + t58 * MDP(15) + (t28 - t32) * MDP(16)) * t42 + t51; (-t32 * t42 - t30) * t62 + (t32 ^ 2 + t33 ^ 2) * MDP(17) + (t47 + t54) * t37 + 0.2e1 * t49 * pkin(2) + t51; -t28 * MDP(14) - t29 * MDP(15) + t48 * t35 + (-t28 * MDP(17) - t53) * pkin(3) + t57; t32 * MDP(14) + t33 * MDP(15) + t48 * pkin(6) + (MDP(17) * t32 - t53) * pkin(3) + t57; MDP(11) + (MDP(17) * pkin(3) + 0.2e1 * MDP(14)) * pkin(3); t50 + t55; t50 + t54; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
