% Calculate joint inertia matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaJ_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRP1_inertiaJ_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:52
% EndTime: 2019-03-08 18:06:52
% DurationCPUTime: 0.03s
% Computational Cost: add. (28->22), mult. (44->27), div. (0->0), fcn. (13->2), ass. (0->10)
t19 = 2 * MDP(7);
t18 = 2 * qJ(3);
t14 = sin(qJ(2));
t17 = t14 * MDP(6);
t16 = pkin(2) * t19 + MDP(4);
t15 = cos(qJ(2));
t13 = t14 * pkin(1);
t11 = t15 * pkin(1) + pkin(2);
t10 = t13 + qJ(3);
t1 = [MDP(1) + MDP(4) + (t10 ^ 2 + t11 ^ 2) * MDP(9) + 0.2e1 * (t15 * MDP(5) - t17) * pkin(1) + t11 * t19 + 0.2e1 * t10 * MDP(8); (t18 + t13) * MDP(8) + (t11 * pkin(2) + t10 * qJ(3)) * MDP(9) + (-t17 + (MDP(5) + MDP(7)) * t15) * pkin(1) + t16; MDP(8) * t18 + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(9) + t16; -t11 * MDP(9) - MDP(7); -pkin(2) * MDP(9) - MDP(7); MDP(9);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
