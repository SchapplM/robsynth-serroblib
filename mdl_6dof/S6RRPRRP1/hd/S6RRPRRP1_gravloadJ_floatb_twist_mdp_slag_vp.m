% Calculate Gravitation load on the joints for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:15
% EndTime: 2019-03-09 11:41:16
% DurationCPUTime: 0.26s
% Computational Cost: add. (280->62), mult. (271->80), div. (0->0), fcn. (230->10), ass. (0->35)
t111 = MDP(19) - MDP(27);
t87 = cos(qJ(2));
t77 = t87 * pkin(2);
t80 = qJ(2) + pkin(10);
t76 = qJ(4) + t80;
t72 = sin(t76);
t73 = cos(t76);
t86 = cos(qJ(5));
t74 = t86 * pkin(5) + pkin(4);
t81 = -qJ(6) - pkin(9);
t94 = -t72 * t81 + t73 * t74;
t110 = t94 + pkin(3) * cos(t80) + t77;
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t70 = g(1) * t88 + g(2) * t85;
t59 = -g(3) * t73 + t70 * t72;
t104 = g(3) * t72;
t83 = sin(qJ(5));
t101 = t85 * t83;
t100 = t85 * t86;
t99 = t88 * t83;
t98 = t88 * t86;
t82 = -qJ(3) - pkin(7);
t95 = pkin(5) * t83 + pkin(8) - t82;
t93 = t111 * (t70 * t73 + t104) + (MDP(25) * t86 - MDP(26) * t83 + MDP(18)) * t59;
t69 = g(1) * t85 - g(2) * t88;
t92 = t72 * t74 + t73 * t81;
t64 = -t73 * t99 + t100;
t62 = t73 * t101 + t98;
t91 = -pkin(1) - t110;
t84 = sin(qJ(2));
t75 = t77 + pkin(1);
t65 = t73 * t98 + t101;
t63 = -t73 * t100 + t99;
t1 = [(-g(1) * (-t85 * t75 - t88 * t82) - g(2) * (t88 * t75 - t85 * t82)) * MDP(12) + (-g(1) * t63 - g(2) * t65) * MDP(25) + (-g(1) * t62 - g(2) * t64) * MDP(26) + ((-g(1) * t95 + g(2) * t91) * t88 + (-g(1) * t91 - g(2) * t95) * t85) * MDP(28) + (MDP(3) - MDP(11)) * t70 + (-t84 * MDP(10) + t73 * MDP(18) + t87 * MDP(9) - t111 * t72 + MDP(2)) * t69; (g(3) * t84 + t70 * t87) * MDP(10) + (-g(3) * t110 + t70 * (pkin(3) * sin(t80) + t84 * pkin(2) + t92)) * MDP(28) + t93 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t87 + t70 * t84); (-MDP(12) - MDP(28)) * t69; (-g(3) * t94 + t70 * t92) * MDP(28) + t93; (g(1) * t65 - g(2) * t63 + t86 * t104) * MDP(26) + (MDP(28) * pkin(5) + MDP(25)) * (-g(1) * t64 + g(2) * t62 + t83 * t104); -t59 * MDP(28);];
taug  = t1;
