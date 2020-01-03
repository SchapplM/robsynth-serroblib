% Calculate joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR11_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (234->73), mult. (416->100), div. (0->0), fcn. (420->6), ass. (0->33)
t59 = sin(pkin(8));
t60 = cos(pkin(8));
t62 = sin(qJ(3));
t79 = cos(qJ(3));
t49 = t62 * t59 - t79 * t60;
t64 = -pkin(3) - pkin(4);
t50 = t79 * t59 + t62 * t60;
t56 = -t60 * pkin(2) - pkin(1);
t67 = t50 * qJ(4) - t56;
t80 = 0.2e1 * t64 * t49 + 0.2e1 * t67;
t78 = pkin(1) * MDP(7);
t77 = pkin(6) + qJ(2);
t75 = t59 * MDP(5);
t74 = t60 * MDP(4);
t61 = sin(qJ(5));
t63 = cos(qJ(5));
t43 = -t63 * t49 + t61 * t50;
t73 = t43 * MDP(24);
t72 = (t61 * qJ(4) - t63 * t64) * MDP(24);
t71 = (t63 * qJ(4) + t61 * t64) * MDP(25);
t70 = MDP(14) - MDP(17);
t53 = t77 * t59;
t54 = t77 * t60;
t45 = t79 * t53 + t62 * t54;
t69 = -pkin(3) * MDP(18) - MDP(15);
t68 = t63 * MDP(24) - t61 * MDP(25);
t46 = -t62 * t53 + t79 * t54;
t40 = -t50 * pkin(7) + t45;
t41 = t49 * pkin(7) + t46;
t44 = t61 * t49 + t63 * t50;
t66 = t44 * MDP(21) - t43 * MDP(22) - (-t63 * t40 + t61 * t41) * MDP(24) - (t61 * t40 + t63 * t41) * MDP(25);
t42 = t49 * pkin(3) - t67;
t1 = [MDP(1) + (t42 ^ 2 + t45 ^ 2 + t46 ^ 2) * MDP(18) + t73 * t80 + (0.2e1 * t74 - 0.2e1 * t75 + t78) * pkin(1) + 0.2e1 * (t56 * MDP(13) + t42 * MDP(15) - MDP(16) * t46) * t49 + (0.2e1 * t56 * MDP(14) + 0.2e1 * MDP(16) * t45 - 0.2e1 * t42 * MDP(17) + MDP(8) * t50 - 0.2e1 * t49 * MDP(9)) * t50 + (MDP(19) * t44 - 0.2e1 * t43 * MDP(20) + MDP(25) * t80) * t44 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t59 ^ 2 + t60 ^ 2) * qJ(2); t42 * MDP(18) - t73 - t44 * MDP(25) - t74 + t75 - t78 + t70 * t50 + (MDP(13) + MDP(15)) * t49; MDP(7) + MDP(18); t50 * MDP(10) - t49 * MDP(11) + (-pkin(3) * t50 - t49 * qJ(4)) * MDP(16) + (MDP(18) * qJ(4) - t70) * t46 + (-MDP(13) + t69) * t45 - t66; 0; MDP(12) + 0.2e1 * pkin(3) * MDP(15) + 0.2e1 * qJ(4) * MDP(17) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(18) + MDP(23) + 0.2e1 * t72 + 0.2e1 * t71; t50 * MDP(16) + t45 * MDP(18); 0; -t68 + t69; MDP(18); t66; 0; -MDP(23) - t71 - t72; t68; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
