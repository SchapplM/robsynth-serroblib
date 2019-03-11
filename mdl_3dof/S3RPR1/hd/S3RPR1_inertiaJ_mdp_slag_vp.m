% Calculate joint inertia matrix for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_inertiaJ_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RPR1_inertiaJ_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.02s
% Computational Cost: add. (21->15), mult. (28->18), div. (0->0), fcn. (12->2), ass. (0->7)
t10 = -pkin(1) - pkin(2);
t8 = sin(qJ(3));
t9 = cos(qJ(3));
t13 = (t8 * qJ(2) - t9 * t10) * MDP(8);
t12 = (t9 * qJ(2) + t8 * t10) * MDP(9);
t11 = t9 * MDP(8) - t8 * MDP(9);
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + MDP(7) + 0.2e1 * t13 + 0.2e1 * t12; -pkin(1) * MDP(6) - MDP(4) - t11; MDP(6); -MDP(7) - t12 - t13; t11; MDP(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
