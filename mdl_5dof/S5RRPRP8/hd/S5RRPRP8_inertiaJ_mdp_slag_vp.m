% Calculate joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP8_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:28
% EndTime: 2019-12-31 20:04:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (195->79), mult. (320->114), div. (0->0), fcn. (266->4), ass. (0->30)
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t48 = -t61 * pkin(2) - t59 * qJ(3) - pkin(1);
t41 = t61 * pkin(3) - t48;
t72 = 0.2e1 * t41;
t71 = pkin(6) - pkin(7);
t56 = t59 ^ 2;
t70 = t61 ^ 2 + t56;
t69 = MDP(23) * pkin(4);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t62 = -pkin(2) - pkin(3);
t46 = qJ(3) * t58 - t60 * t62;
t68 = MDP(20) * t46;
t47 = qJ(3) * t60 + t58 * t62;
t67 = MDP(21) * t47;
t66 = t58 * MDP(21);
t49 = t71 * t59;
t50 = t71 * t61;
t39 = -t60 * t49 + t50 * t58;
t65 = -pkin(2) * MDP(14) - MDP(11);
t40 = t49 * t58 + t50 * t60;
t42 = t58 * t59 + t60 * t61;
t43 = -t61 * t58 + t59 * t60;
t64 = t43 * MDP(17) - t42 * MDP(18) - t39 * MDP(20) - t40 * MDP(21);
t45 = -pkin(4) - t46;
t38 = pkin(4) * t42 + t41;
t37 = -qJ(5) * t42 + t40;
t36 = -qJ(5) * t43 - t39;
t1 = [MDP(1) + t56 * MDP(4) + (t70 * pkin(6) ^ 2 + t48 ^ 2) * MDP(14) + t42 * MDP(20) * t72 + (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) * MDP(23) + (MDP(15) * t43 - 0.2e1 * t42 * MDP(16) + MDP(21) * t72) * t43 + 0.2e1 * (-t36 * t43 - t37 * t42) * MDP(22) + 0.2e1 * t70 * MDP(12) * pkin(6) + 0.2e1 * (-t48 * MDP(11) + pkin(1) * MDP(9)) * t61 + 0.2e1 * (-MDP(10) * pkin(1) - MDP(13) * t48 + MDP(5) * t61) * t59; t59 * MDP(6) + t61 * MDP(7) + (-pkin(2) * t59 + qJ(3) * t61) * MDP(12) + (-t42 * t47 - t43 * t45) * MDP(22) + (t36 * t45 + t37 * t47) * MDP(23) + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t61 + (-MDP(9) + t65) * t59) * pkin(6) - t64; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + MDP(19) + (t45 ^ 2 + t47 ^ 2) * MDP(23) + 0.2e1 * t68 + 0.2e1 * t67; (-t42 * t58 - t43 * t60) * MDP(22) + (t36 * t60 + t37 * t58) * MDP(23) + (MDP(14) * pkin(6) + MDP(12)) * t59; -t60 * MDP(20) + t66 + (t45 * t60 + t47 * t58) * MDP(23) + t65; MDP(14) + (t58 ^ 2 + t60 ^ 2) * MDP(23); (-MDP(22) * t43 + MDP(23) * t36) * pkin(4) + t64; t45 * t69 - MDP(19) - t67 - t68; -t66 + (MDP(20) + t69) * t60; pkin(4) ^ 2 * MDP(23) + MDP(19); t38 * MDP(23); 0; 0; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
