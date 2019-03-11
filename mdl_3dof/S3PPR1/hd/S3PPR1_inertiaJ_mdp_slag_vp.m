% Calculate joint inertia matrix for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% MDP [5x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3PPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_inertiaJ_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_inertiaJ_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [5 1]), ...
  'S3PPR1_inertiaJ_mdp_slag_vp: MDP has to be [5x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:02:00
% EndTime: 2019-03-08 18:02:00
% DurationCPUTime: 0.01s
% Computational Cost: add. (3->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
t4 = cos(qJ(3));
t3 = sin(qJ(3));
t1 = [MDP(1) + MDP(2); 0; MDP(2); -t3 * MDP(4) - t4 * MDP(5); t4 * MDP(4) - t3 * MDP(5); MDP(3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
