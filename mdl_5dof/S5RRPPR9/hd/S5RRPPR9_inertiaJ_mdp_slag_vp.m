% Calculate joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR9_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:48
% EndTime: 2019-12-31 19:41:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (163->84), mult. (269->107), div. (0->0), fcn. (177->4), ass. (0->33)
t80 = pkin(6) - qJ(4);
t64 = -pkin(2) - pkin(3);
t63 = cos(qJ(2));
t48 = t80 * t63;
t79 = t48 * t63;
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t78 = t60 * t62;
t77 = pkin(6) ^ 2 * MDP(14);
t76 = MDP(12) - MDP(17);
t75 = 0.2e1 * t63;
t61 = sin(qJ(2));
t46 = -t63 * pkin(2) - t61 * qJ(3) - pkin(1);
t74 = MDP(20) * t78;
t45 = t63 * pkin(3) - t46;
t73 = -pkin(2) * MDP(14) - MDP(11);
t72 = t64 * MDP(18) + MDP(16);
t71 = -MDP(21) * t62 + MDP(22) * t60;
t70 = t62 * MDP(24) - t60 * MDP(25);
t69 = -t60 * MDP(24) - t62 * MDP(25);
t68 = MDP(15) + t70;
t67 = -t60 * MDP(21) - t62 * MDP(22) + t69 * (-pkin(7) + t64);
t65 = qJ(3) ^ 2;
t59 = qJ(3) + pkin(4);
t58 = t63 ^ 2;
t57 = t62 ^ 2;
t56 = t61 ^ 2;
t55 = t60 ^ 2;
t47 = t80 * t61;
t44 = t61 * pkin(4) + t63 * pkin(7) + t45;
t43 = t60 * t44 + t62 * t47;
t42 = t62 * t44 - t60 * t47;
t1 = [MDP(1) + t46 ^ 2 * MDP(14) + (t45 ^ 2 + t47 ^ 2 + t48 ^ 2) * MDP(18) + (t57 * MDP(19) - 0.2e1 * t74 + t77) * t58 + (MDP(23) + MDP(4) + t77) * t56 + (-t46 * MDP(11) - t45 * MDP(16) + pkin(1) * MDP(9)) * t75 + (-0.2e1 * pkin(1) * MDP(10) - 0.2e1 * t46 * MDP(13) + 0.2e1 * t45 * MDP(15) + (MDP(5) + t71) * t75) * t61 + 0.2e1 * (-t47 * t61 - t79) * MDP(17) + 0.2e1 * (t42 * t61 - t60 * t79) * MDP(24) + 0.2e1 * (-t43 * t61 - t62 * t79) * MDP(25) + 0.2e1 * (t56 + t58) * MDP(12) * pkin(6); t72 * t47 + (qJ(3) * MDP(18) + t68) * t48 + (-pkin(2) * MDP(12) - t64 * MDP(17) + MDP(6) + t67) * t61 + (MDP(7) + MDP(19) * t78 + (-t55 + t57) * MDP(20) + t69 * t59 + t76 * qJ(3)) * t63 + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t63 + (-MDP(9) + t73) * t61) * pkin(6); MDP(8) + 0.2e1 * pkin(2) * MDP(11) + (pkin(2) ^ 2 + t65) * MDP(14) + 0.2e1 * t64 * MDP(16) + (t64 ^ 2 + t65) * MDP(18) + t55 * MDP(19) + 0.2e1 * t74 + 0.2e1 * t70 * t59 + 0.2e1 * (MDP(13) + MDP(15)) * qJ(3); t47 * MDP(18) + (MDP(14) * pkin(6) + t69 + t76) * t61; t72 + t73; MDP(14) + MDP(18); -t63 * MDP(16) + t45 * MDP(18) + t61 * t68; 0; 0; MDP(18); t61 * MDP(23) + t42 * MDP(24) - t43 * MDP(25) + t63 * t71; t67; t69; t70; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
