% Calculate joint inertia matrix for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% MDP [7x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3PRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_inertiaJ_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [7 1]), ...
  'S3PRR1_inertiaJ_mdp_slag_vp: MDP has to be [7x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:00
% EndTime: 2019-03-08 18:04:01
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->9), mult. (24->12), div. (0->0), fcn. (22->4), ass. (0->7)
t11 = sin(qJ(3));
t12 = sin(qJ(2));
t13 = cos(qJ(3));
t14 = cos(qJ(2));
t16 = (-t11 * t12 + t13 * t14) * MDP(6) + (-t11 * t14 - t13 * t12) * MDP(7);
t15 = (MDP(6) * t13 - MDP(7) * t11) * pkin(2);
t1 = [MDP(1); t14 * MDP(3) - t12 * MDP(4) + t16; MDP(2) + MDP(5) + 0.2e1 * t15; t16; MDP(5) + t15; MDP(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
