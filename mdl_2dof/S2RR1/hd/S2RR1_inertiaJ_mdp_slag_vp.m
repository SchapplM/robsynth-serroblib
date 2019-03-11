% Calculate joint inertia matrix for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_inertiaJ_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_inertiaJ_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S2RR1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:03
% EndTime: 2019-03-08 18:00:03
% DurationCPUTime: 0.02s
% Computational Cost: add. (5->5), mult. (11->9), div. (0->0), fcn. (7->2), ass. (0->3)
t4 = cos(qJ(2));
t3 = sin(qJ(2));
t1 = [MDP(1) + (MDP(4) * t3 + 0.2e1 * MDP(5) * t4) * t3; -t3 * MDP(6) - t4 * MDP(7) + (MDP(10) * t4 + MDP(9) * t3) * pkin(1); MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t1(1) t1(2); t1(2) t1(3);];
Mq  = res;
