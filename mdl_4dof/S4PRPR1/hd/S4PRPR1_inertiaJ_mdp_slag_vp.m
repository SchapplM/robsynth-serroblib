% Calculate joint inertia matrix for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPR1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:03
% EndTime: 2019-03-08 18:21:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->16), mult. (28->18), div. (0->0), fcn. (12->2), ass. (0->7)
t10 = -pkin(2) - pkin(3);
t8 = sin(qJ(4));
t9 = cos(qJ(4));
t13 = (t8 * qJ(3) - t9 * t10) * MDP(9);
t12 = (t9 * qJ(3) + t8 * t10) * MDP(10);
t11 = -t8 * MDP(10) + t9 * MDP(9);
t1 = [MDP(1) + MDP(7); 0; MDP(2) + (2 * pkin(2) * MDP(5)) + 0.2e1 * qJ(3) * MDP(6) + ((pkin(2) ^ 2) + qJ(3) ^ 2) * MDP(7) + MDP(8) + 0.2e1 * t12 + 0.2e1 * t13; 0; -pkin(2) * MDP(7) - MDP(5) - t11; MDP(7); 0; -MDP(8) - t12 - t13; t11; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
