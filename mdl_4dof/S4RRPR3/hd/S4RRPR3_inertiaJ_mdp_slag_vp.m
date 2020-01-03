% Calculate joint inertia matrix for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RRPR3_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->28), mult. (118->40), div. (0->0), fcn. (81->6), ass. (0->19)
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t47 = t41 * MDP(13) - t39 * MDP(14);
t40 = sin(qJ(2));
t53 = pkin(1) * t40;
t42 = cos(qJ(2));
t33 = t42 * pkin(1) + pkin(2);
t37 = sin(pkin(7));
t38 = cos(pkin(7));
t27 = t37 * t33 + t38 * t53;
t51 = t39 * MDP(10) + t41 * MDP(11);
t48 = MDP(4) + (MDP(8) * t39 + 0.2e1 * MDP(9) * t41) * t39;
t46 = -MDP(13) * t39 - MDP(14) * t41;
t26 = t38 * t33 - t37 * t53;
t45 = (t42 * MDP(5) - t40 * MDP(6)) * pkin(1);
t44 = 0.2e1 * t47;
t32 = -t38 * pkin(2) - pkin(3);
t24 = -pkin(3) - t26;
t1 = [MDP(1) + (t26 ^ 2 + t27 ^ 2) * MDP(7) - t24 * t44 + 0.2e1 * t45 + t48; (t26 * t38 + t27 * t37) * MDP(7) * pkin(2) + t45 + t48 - t47 * (t24 + t32); (t37 ^ 2 + t38 ^ 2) * MDP(7) * pkin(2) ^ 2 - t32 * t44 + t48; 0; 0; MDP(7); t46 * (pkin(6) + t27) + t51; t46 * (t37 * pkin(2) + pkin(6)) + t51; t47; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
