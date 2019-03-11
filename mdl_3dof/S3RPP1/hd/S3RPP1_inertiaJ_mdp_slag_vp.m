% Calculate joint inertia matrix for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_inertiaJ_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_inertiaJ_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPP1_inertiaJ_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:02
% EndTime: 2019-03-08 18:05:02
% DurationCPUTime: 0.02s
% Computational Cost: add. (16->14), mult. (18->14), div. (0->0), fcn. (0->0), ass. (0->3)
t7 = (qJ(2) ^ 2);
t5 = (pkin(1) + qJ(3));
t1 = [MDP(1) - 2 * pkin(1) * MDP(4) + (pkin(1) ^ 2 + t7) * MDP(6) + 2 * t5 * MDP(8) + (t5 ^ 2 + t7) * MDP(9) + 2 * (MDP(5) + MDP(7)) * qJ(2); -pkin(1) * MDP(6) - t5 * MDP(9) + MDP(4) - MDP(8); MDP(6) + MDP(9); qJ(2) * MDP(9) + MDP(7); 0; MDP(9);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
