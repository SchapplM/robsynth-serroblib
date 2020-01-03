% Calculate joint inertia matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR4_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:24
% EndTime: 2019-12-31 17:32:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (60->31), mult. (128->51), div. (0->0), fcn. (111->6), ass. (0->22)
t33 = cos(pkin(8));
t42 = t33 * MDP(6);
t32 = sin(pkin(8));
t43 = t32 * MDP(7);
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t23 = t34 * t32 - t36 * t33;
t21 = t23 * MDP(15);
t24 = t36 * t32 + t34 * t33;
t45 = t24 * MDP(16) + t21;
t47 = pkin(3) * MDP(9);
t49 = -t42 + t43 - t47 + t45;
t48 = -0.2e1 * t33 * pkin(4) - (2 * pkin(3));
t46 = pkin(6) + qJ(4);
t44 = t32 ^ 2 + t33 ^ 2;
t41 = t44 * MDP(9);
t40 = t44 * qJ(4);
t37 = cos(qJ(3));
t35 = sin(qJ(3));
t26 = t46 * t33;
t25 = t46 * t32;
t1 = [MDP(1) + MDP(2) + t41; 0; MDP(2) + (t44 * t35 ^ 2 + t37 ^ 2) * MDP(9); 0; (t44 * MDP(8) + MDP(9) * t40 - MDP(5)) * t35 + (MDP(4) - t49) * t37; t21 * t48 + MDP(3) + qJ(4) ^ 2 * t41 + 0.2e1 * MDP(8) * t40 + (0.2e1 * t42 - 0.2e1 * t43 + t47) * pkin(3) + (MDP(10) * t24 - 0.2e1 * t23 * MDP(11) + MDP(16) * t48) * t24; 0; -t37 * MDP(9); t49; MDP(9); t45; (-t24 * MDP(15) + t23 * MDP(16)) * t35; t24 * MDP(12) - t23 * MDP(13) + (-t36 * t25 - t34 * t26) * MDP(15) + (t34 * t25 - t36 * t26) * MDP(16); 0; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
