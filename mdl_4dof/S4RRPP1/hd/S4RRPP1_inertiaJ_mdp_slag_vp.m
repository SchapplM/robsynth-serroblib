% Calculate joint inertia matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPP1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:10
% EndTime: 2019-03-08 18:33:10
% DurationCPUTime: 0.05s
% Computational Cost: add. (70->38), mult. (115->50), div. (0->0), fcn. (69->4), ass. (0->17)
t41 = 2 * MDP(9);
t36 = sin(qJ(2));
t40 = pkin(1) * t36;
t37 = cos(qJ(2));
t32 = t37 * pkin(1) + pkin(2);
t34 = sin(pkin(6));
t35 = cos(pkin(6));
t25 = t34 * t32 + t35 * t40;
t39 = t37 * MDP(5);
t27 = t35 * t32;
t24 = -t34 * t40 + t27;
t33 = t34 * pkin(2);
t30 = t35 * pkin(2) + pkin(3);
t29 = t33 + qJ(4);
t23 = -pkin(3) - t24;
t22 = qJ(4) + t25;
t1 = [MDP(1) + MDP(4) + (t24 ^ 2 + t25 ^ 2) * MDP(7) + (t22 ^ 2 + t23 ^ 2) * MDP(10) + 0.2e1 * (-t36 * MDP(6) + t39) * pkin(1) - 0.2e1 * t23 * MDP(8) + t22 * t41; MDP(4) + (0.2e1 * pkin(3) + t27) * MDP(8) + (t33 + 0.2e1 * qJ(4) + t25) * MDP(9) + (t22 * t29 - t23 * t30) * MDP(10) + ((t24 * t35 + t25 * t34) * MDP(7) + t35 * MDP(8)) * pkin(2) + (t39 + (-MDP(8) * t34 - MDP(6)) * t36) * pkin(1); MDP(4) + (t29 ^ 2 + t30 ^ 2) * MDP(10) + (t34 ^ 2 + t35 ^ 2) * MDP(7) * pkin(2) ^ 2 + 0.2e1 * t30 * MDP(8) + t29 * t41; 0; 0; MDP(7) + MDP(10); t23 * MDP(10) - MDP(8); -t30 * MDP(10) - MDP(8); 0; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
