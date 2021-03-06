% Calculate joint inertia matrix for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR4_inertiaJ_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:40
% DurationCPUTime: 0.03s
% Computational Cost: add. (20->14), mult. (44->25), div. (0->0), fcn. (41->6), ass. (0->11)
t16 = cos(qJ(4));
t20 = t16 * MDP(11);
t14 = sin(qJ(4));
t19 = -t14 * MDP(12) + t20;
t18 = -MDP(11) * t14 - MDP(12) * t16;
t17 = cos(qJ(3));
t15 = sin(qJ(3));
t13 = cos(pkin(7));
t12 = sin(pkin(7));
t11 = t17 * t12 + t15 * t13;
t1 = [MDP(1) + (t12 ^ 2 + t13 ^ 2) * MDP(2); 0; MDP(2); -t11 * MDP(5) + (-MDP(4) - t19) * (t15 * t12 - t17 * t13); 0; 0.2e1 * pkin(3) * t20 + MDP(3) + (-0.2e1 * MDP(12) * pkin(3) + MDP(6) * t14 + 0.2e1 * MDP(7) * t16) * t14; t18 * t11; t19; t14 * MDP(8) + t16 * MDP(9) + t18 * pkin(5); MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
