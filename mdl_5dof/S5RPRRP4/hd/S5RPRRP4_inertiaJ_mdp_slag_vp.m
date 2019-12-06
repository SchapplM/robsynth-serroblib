% Calculate joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:54
% EndTime: 2019-12-05 18:06:56
% DurationCPUTime: 0.29s
% Computational Cost: add. (310->81), mult. (625->120), div. (0->0), fcn. (589->6), ass. (0->39)
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t65 = -t73 * pkin(2) - t72 * pkin(6) - pkin(1);
t77 = cos(qJ(3));
t61 = t77 * t65;
t75 = sin(qJ(3));
t91 = qJ(2) * t75;
t99 = pkin(7) * t72;
t50 = -t77 * t99 + t61 + (-pkin(3) - t91) * t73;
t84 = t77 * t73 * qJ(2);
t52 = t84 + (t65 - t99) * t75;
t74 = sin(qJ(4));
t76 = cos(qJ(4));
t47 = t76 * t50 - t74 * t52;
t48 = t74 * t50 + t76 * t52;
t64 = t74 * t77 + t76 * t75;
t57 = t64 * t72;
t63 = -t74 * t75 + t76 * t77;
t58 = t63 * t72;
t106 = t58 * MDP(17) - t57 * MDP(18) + t47 * MDP(20) - t48 * MDP(21);
t105 = (t77 * MDP(10) - t75 * MDP(11)) * t72 + (-t73 * t91 + t61) * MDP(13) - (t75 * t65 + t84) * MDP(14) + t106;
t86 = t76 * MDP(20);
t104 = (-t74 * MDP(21) + t86) * pkin(3);
t94 = t63 * MDP(20) - t64 * MDP(21);
t101 = -t77 * MDP(13) + t75 * MDP(14) - t94;
t100 = pkin(3) * t74;
t98 = pkin(1) * MDP(7);
t93 = (pkin(3) * t75 + qJ(2)) * t72;
t92 = MDP(23) * pkin(4);
t90 = t72 * MDP(5);
t89 = qJ(2) ^ 2 * MDP(7);
t85 = MDP(12) + MDP(19);
t71 = t73 ^ 2;
t70 = t72 ^ 2;
t68 = t76 * pkin(3) + pkin(4);
t51 = t57 * pkin(4) + t93;
t46 = -t57 * qJ(5) + t48;
t45 = -t73 * pkin(4) - t58 * qJ(5) + t47;
t1 = [MDP(1) + (t45 ^ 2 + t46 ^ 2 + t51 ^ 2) * MDP(23) + (-0.2e1 * t90 + t98) * pkin(1) + (t85 + t89) * t71 + (t89 + (MDP(8) * t77 - 0.2e1 * MDP(9) * t75) * t77) * t70 + (MDP(15) * t58 - 0.2e1 * t57 * MDP(16)) * t58 + 0.2e1 * (-t45 * t58 - t46 * t57) * MDP(22) + 0.2e1 * (t57 * MDP(20) + t58 * MDP(21)) * t93 + 0.2e1 * (t71 * MDP(6) + (t75 * MDP(13) + t77 * MDP(14) + MDP(6)) * t70) * qJ(2) + 0.2e1 * (pkin(1) * MDP(4) - t105) * t73; t90 - t98 + (-t64 * t57 - t63 * t58) * MDP(22) + (t45 * t63 + t46 * t64) * MDP(23) + (-MDP(4) + t101) * t73; MDP(7) + (t63 ^ 2 + t64 ^ 2) * MDP(23); (-t57 * t100 - t68 * t58) * MDP(22) + (t46 * t100 + t45 * t68) * MDP(23) + (-t85 - t104) * t73 + t105; (t64 * t100 + t63 * t68) * MDP(23) - t101; t68 ^ 2 * MDP(23) + (0.2e1 * t86 + (MDP(23) * t100 - 0.2e1 * MDP(21)) * t74) * pkin(3) + t85; -t73 * MDP(19) + (-t58 * MDP(22) + t45 * MDP(23)) * pkin(4) + t106; t63 * t92 + t94; t68 * t92 + MDP(19) + t104; MDP(23) * pkin(4) ^ 2 + MDP(19); t51 * MDP(23); 0; 0; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
