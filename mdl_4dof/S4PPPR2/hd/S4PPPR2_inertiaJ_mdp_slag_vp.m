% Calculate joint inertia matrix for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPPR2_inertiaJ_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:17
% EndTime: 2019-03-08 18:10:17
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->8), mult. (15->12), div. (0->0), fcn. (15->4), ass. (0->6)
t11 = MDP(2) + MDP(3);
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t8 = cos(pkin(5));
t7 = sin(pkin(5));
t1 = [MDP(1) + t11 * (t7 ^ 2 + t8 ^ 2); 0; t11; -t8 * MDP(3); 0; MDP(3); (-t8 * t10 - t7 * t9) * MDP(5) + (-t7 * t10 + t8 * t9) * MDP(6); 0; t10 * MDP(5) - t9 * MDP(6); MDP(4);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
