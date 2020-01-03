% Calculate joint inertia matrix for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR5_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:31
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->25), mult. (48->26), div. (0->0), fcn. (30->4), ass. (0->12)
t18 = (MDP(8) * pkin(3));
t7 = sin(qJ(5));
t17 = t7 * MDP(14);
t9 = cos(qJ(5));
t16 = t9 * MDP(15);
t15 = qJ(4) * MDP(8);
t14 = MDP(6) - t18;
t13 = t9 * MDP(14) - t7 * MDP(15);
t12 = t16 + t17;
t10 = cos(qJ(3));
t8 = sin(qJ(3));
t1 = [MDP(1) + MDP(2) + MDP(8); 0; MDP(2) + (t10 ^ 2 + t8 ^ 2) * MDP(8); 0; (MDP(4) - t14) * t10 + (-MDP(5) + MDP(7) + t12 + t15) * t8; MDP(3) + (-0.2e1 * t7 * MDP(10) + MDP(9) * t9) * t9 + ((-2 * MDP(6) + t18) * pkin(3)) + (0.2e1 * MDP(7) + t15 + 0.2e1 * t16 + 0.2e1 * t17) * qJ(4); 0; -t10 * MDP(8); t14; MDP(8); t12; -t13 * t10; t9 * MDP(11) - t7 * MDP(12) + t13 * (-pkin(3) - pkin(6)); t13; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
