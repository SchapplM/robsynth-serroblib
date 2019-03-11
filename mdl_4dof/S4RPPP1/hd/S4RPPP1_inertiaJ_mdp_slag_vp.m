% Calculate joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPP1_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:26
% EndTime: 2019-03-08 18:26:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (113->56), mult. (268->76), div. (0->0), fcn. (220->4), ass. (0->23)
t40 = cos(pkin(4));
t48 = pkin(1) * t40;
t37 = sin(pkin(6));
t38 = sin(pkin(4));
t47 = t37 * t38;
t39 = cos(pkin(6));
t46 = t38 * t39;
t45 = qJ(2) * t38;
t32 = t37 * t48 + t39 * t45;
t44 = -MDP(14) + MDP(9);
t43 = MDP(11) + MDP(15);
t42 = -pkin(1) * t39 - pkin(2);
t41 = -qJ(3) * t37 - pkin(1);
t28 = -t40 * qJ(3) - t32;
t36 = t38 ^ 2;
t33 = t37 * t45;
t31 = t39 * t48 - t33;
t30 = (-pkin(2) * t39 + t41) * t38;
t29 = t42 * t40 + t33;
t27 = ((-pkin(2) - qJ(4)) * t39 + t41) * t38;
t26 = pkin(3) * t46 - t28;
t25 = pkin(3) * t47 + t33 + (-qJ(4) + t42) * t40;
t1 = [MDP(1) + (t36 * pkin(1) ^ 2 + t31 ^ 2 + t32 ^ 2) * MDP(7) + (t28 ^ 2 + t29 ^ 2 + t30 ^ 2) * MDP(11) + (t25 ^ 2 + t26 ^ 2 + t27 ^ 2) * MDP(15) + 0.2e1 * (MDP(4) * t39 - MDP(5) * t37) * t36 * pkin(1) + 0.2e1 * (-t28 * MDP(10) + t26 * MDP(13) - t25 * MDP(14) + t31 * MDP(4) - t32 * MDP(5) + t29 * MDP(9)) * t40 + 0.2e1 * ((-t31 * t37 + t32 * t39) * MDP(6) + (-t28 * t39 + t29 * t37) * MDP(8) + (t25 * t37 + t26 * t39) * MDP(12) + (-t37 * MDP(10) + t39 * MDP(9)) * t30 + (-t37 * MDP(13) - t39 * MDP(14)) * t27) * t38; t30 * MDP(11) + t27 * MDP(15) + (-MDP(7) * pkin(1) + (-MDP(4) + t44) * t39 + (-MDP(10) - MDP(13) + MDP(5)) * t37) * t38; MDP(7) + t43; t29 * MDP(11) + t25 * MDP(15) + t44 * t40 + (MDP(12) + MDP(8)) * t47; 0; t43; MDP(12) * t46 + t40 * MDP(13) + t26 * MDP(15); 0; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
