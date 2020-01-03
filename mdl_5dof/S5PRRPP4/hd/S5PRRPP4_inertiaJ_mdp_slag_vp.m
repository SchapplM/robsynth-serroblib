% Calculate joint inertia matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:13
% EndTime: 2019-12-31 17:41:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (100->62), mult. (161->76), div. (0->0), fcn. (90->2), ass. (0->20)
t49 = pkin(6) - qJ(5);
t35 = sin(qJ(3));
t36 = cos(qJ(3));
t48 = t36 * MDP(16) + t35 * MDP(17);
t31 = t35 * qJ(4);
t47 = t36 * pkin(3) + t31;
t33 = t35 ^ 2;
t46 = t36 ^ 2 + t33;
t45 = pkin(3) * MDP(15);
t44 = t36 * qJ(4);
t43 = -MDP(11) + MDP(14);
t42 = MDP(15) + MDP(19);
t25 = -pkin(2) - t47;
t37 = pkin(3) + pkin(4);
t41 = t37 * MDP(19) + MDP(12);
t39 = qJ(4) ^ 2;
t27 = t49 * t36;
t26 = t49 * t35;
t24 = t36 * pkin(4) - t25;
t1 = [t42 * t46 + MDP(1); (-t36 * t26 + t35 * t27) * MDP(19); MDP(2) + t33 * MDP(5) + (t46 * pkin(6) ^ 2 + t25 ^ 2) * MDP(15) + (t24 ^ 2 + t26 ^ 2 + t27 ^ 2) * MDP(19) + 0.2e1 * t46 * MDP(13) * pkin(6) + 0.2e1 * (pkin(2) * MDP(10) - t25 * MDP(12) + t24 * MDP(16) - MDP(18) * t27) * t36 + 0.2e1 * (-pkin(2) * MDP(11) - t25 * MDP(14) + t24 * MDP(17) - MDP(18) * t26 + t36 * MDP(6)) * t35; t47 * MDP(15) + t31 * MDP(19) + (MDP(10) + t41) * t36 + t43 * t35 + t48; t35 * MDP(7) + t36 * MDP(8) + (-t35 * pkin(3) + t44) * MDP(13) - t26 * MDP(16) + t27 * MDP(17) + (t37 * t35 - t44) * MDP(18) + (t27 * qJ(4) - t26 * t37) * MDP(19) + ((qJ(4) * MDP(15) + t43) * t36 + (-MDP(10) - MDP(12) - t45) * t35) * pkin(6); MDP(9) + 0.2e1 * pkin(3) * MDP(12) + (pkin(3) ^ 2 + t39) * MDP(15) + 0.2e1 * t37 * MDP(16) + (t37 ^ 2 + t39) * MDP(19) + 0.2e1 * (MDP(14) + MDP(17)) * qJ(4); -t42 * t36; t26 * MDP(19) + (MDP(15) * pkin(6) + MDP(13) - MDP(18)) * t35; -MDP(16) - t41 - t45; t42; 0; t24 * MDP(19) + t48; 0; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
