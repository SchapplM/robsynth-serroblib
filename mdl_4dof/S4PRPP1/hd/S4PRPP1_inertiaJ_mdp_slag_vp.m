% Calculate joint inertia matrix for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:55
% EndTime: 2019-03-08 18:17:55
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->15), mult. (18->14), div. (0->0), fcn. (0->0), ass. (0->4)
t8 = MDP(7) + MDP(10);
t7 = (qJ(3) ^ 2);
t5 = (pkin(2) + qJ(4));
t1 = [MDP(1) + t8; 0; MDP(2) - 2 * pkin(2) * MDP(5) + (pkin(2) ^ 2 + t7) * MDP(7) + 2 * t5 * MDP(9) + (t5 ^ 2 + t7) * MDP(10) + 2 * (MDP(6) + MDP(8)) * qJ(3); 0; -t5 * MDP(10) - pkin(2) * MDP(7) + MDP(5) - MDP(9); t8; 0; qJ(3) * MDP(10) + MDP(8); 0; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
