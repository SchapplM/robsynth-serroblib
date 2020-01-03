% Calculate Gravitation load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:17
% EndTime: 2019-12-31 21:17:18
% DurationCPUTime: 0.29s
% Computational Cost: add. (224->62), mult. (265->90), div. (0->0), fcn. (236->10), ass. (0->39)
t80 = cos(qJ(2));
t75 = qJ(2) + qJ(3);
t71 = sin(t75);
t72 = cos(t75);
t90 = t72 * pkin(3) + t71 * qJ(4);
t103 = t80 * pkin(2) + t90;
t102 = MDP(17) - MDP(20);
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t87 = g(1) * t81 + g(2) * t79;
t57 = -g(3) * t72 + t87 * t71;
t101 = pkin(3) * t71;
t100 = g(3) * t71;
t74 = pkin(9) + qJ(5);
t69 = sin(t74);
t98 = t79 * t69;
t70 = cos(t74);
t97 = t79 * t70;
t76 = sin(pkin(9));
t96 = t79 * t76;
t77 = cos(pkin(9));
t95 = t79 * t77;
t94 = t81 * t69;
t93 = t81 * t70;
t92 = t81 * t76;
t91 = t81 * t77;
t89 = qJ(4) * t72;
t78 = sin(qJ(2));
t88 = -pkin(2) * t78 - t101;
t85 = pkin(1) + t103;
t84 = t102 * (t87 * t72 + t100) + (MDP(18) * t77 - MDP(19) * t76 + MDP(27) * t70 - MDP(28) * t69 + MDP(16)) * t57;
t82 = -pkin(7) - pkin(6);
t65 = t81 * t89;
t64 = t79 * t89;
t62 = t72 * t93 + t98;
t61 = -t72 * t94 + t97;
t60 = -t72 * t97 + t94;
t59 = t72 * t98 + t93;
t1 = [t87 * MDP(3) + (-g(1) * (-t72 * t95 + t92) - g(2) * (t72 * t91 + t96)) * MDP(18) + (-g(1) * (t72 * t96 + t91) - g(2) * (-t72 * t92 + t95)) * MDP(19) + ((g(1) * t82 - g(2) * t85) * t81 + (g(1) * t85 + g(2) * t82) * t79) * MDP(21) + (-g(1) * t60 - g(2) * t62) * MDP(27) + (-g(1) * t59 - g(2) * t61) * MDP(28) + (-t78 * MDP(10) + t72 * MDP(16) + t80 * MDP(9) - t102 * t71 + MDP(2)) * (g(1) * t79 - g(2) * t81); (-g(3) * t80 + t87 * t78) * MDP(9) + (g(3) * t78 + t87 * t80) * MDP(10) + (-g(1) * (t88 * t81 + t65) - g(2) * (t88 * t79 + t64) - g(3) * t103) * MDP(21) + t84; (-g(1) * (-t81 * t101 + t65) - g(2) * (-t79 * t101 + t64) - g(3) * t90) * MDP(21) + t84; -t57 * MDP(21); (-g(1) * t61 + g(2) * t59 + t69 * t100) * MDP(27) + (g(1) * t62 - g(2) * t60 + t70 * t100) * MDP(28);];
taug = t1;
