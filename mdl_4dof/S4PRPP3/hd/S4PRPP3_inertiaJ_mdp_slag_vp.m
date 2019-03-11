% Calculate joint inertia matrix for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRPP3_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:54
% EndTime: 2019-03-08 18:19:54
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->20), mult. (37->20), div. (0->0), fcn. (16->2), ass. (0->8)
t18 = MDP(6) + MDP(9);
t17 = MDP(7) + MDP(10);
t13 = (pkin(2) + pkin(3));
t16 = -t13 * MDP(10) - pkin(2) * MDP(7) - MDP(5) - MDP(8);
t15 = (qJ(3) ^ 2);
t12 = cos(qJ(2));
t11 = sin(qJ(2));
t1 = [MDP(1) + t17 * (t11 ^ 2 + t12 ^ 2); (MDP(3) - t16) * t12 + (t17 * qJ(3) - MDP(4) + t18) * t11; MDP(2) + (2 * pkin(2) * MDP(5)) + ((pkin(2) ^ 2 + t15) * MDP(7)) + (2 * t13 * MDP(8)) + ((t13 ^ 2 + t15) * MDP(10)) + 0.2e1 * t18 * qJ(3); -t17 * t12; t16; t17; 0; 0; 0; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
