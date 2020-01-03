% Calculate Gravitation load on the joints for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:39
% EndTime: 2019-12-31 20:38:40
% DurationCPUTime: 0.47s
% Computational Cost: add. (259->88), mult. (492->152), div. (0->0), fcn. (574->12), ass. (0->39)
t72 = sin(qJ(2));
t73 = sin(qJ(1));
t75 = cos(qJ(2));
t76 = cos(qJ(1));
t84 = cos(pkin(5));
t82 = t76 * t84;
t59 = t72 * t82 + t73 * t75;
t67 = pkin(10) + qJ(4);
t65 = sin(t67);
t66 = cos(t67);
t69 = sin(pkin(5));
t87 = t69 * t76;
t52 = t59 * t66 - t65 * t87;
t58 = t73 * t72 - t75 * t82;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t97 = t52 * t71 - t58 * t74;
t96 = t52 * t74 + t58 * t71;
t95 = MDP(10) - MDP(13);
t94 = g(3) * t69;
t91 = t66 * t71;
t90 = t66 * t74;
t89 = t69 * t72;
t88 = t69 * t73;
t86 = t71 * t75;
t85 = t74 * t75;
t83 = t73 * t84;
t80 = t59 * t65 + t66 * t87;
t60 = t76 * t72 + t75 * t83;
t78 = -g(1) * t60 - g(2) * t58 + t75 * t94;
t70 = cos(pkin(10));
t68 = sin(pkin(10));
t61 = -t72 * t83 + t76 * t75;
t57 = t84 * t65 + t66 * t89;
t55 = t61 * t66 + t65 * t88;
t54 = -t61 * t65 + t66 * t88;
t50 = t55 * t74 + t60 * t71;
t49 = -t55 * t71 + t60 * t74;
t1 = [(g(1) * t73 - g(2) * t76) * MDP(2) + (g(1) * t76 + g(2) * t73) * MDP(3) + (g(1) * t59 - g(2) * t61) * MDP(9) + (-g(1) * (-t59 * t70 + t68 * t87) - g(2) * (t61 * t70 + t68 * t88)) * MDP(11) + (-g(1) * (t59 * t68 + t70 * t87) - g(2) * (-t61 * t68 + t70 * t88)) * MDP(12) + (-g(1) * (-t73 * pkin(1) - t59 * pkin(2) + pkin(7) * t87 - t58 * qJ(3)) - g(2) * (t76 * pkin(1) + t61 * pkin(2) + pkin(7) * t88 + t60 * qJ(3))) * MDP(14) + (g(1) * t52 - g(2) * t55) * MDP(20) + (-g(1) * t80 - g(2) * t54) * MDP(21) + (g(1) * t96 - g(2) * t50) * MDP(27) + (-g(1) * t97 - g(2) * t49) * MDP(28) - t95 * (g(1) * t58 - g(2) * t60); (-g(1) * (-t60 * pkin(2) + t61 * qJ(3)) - g(2) * (-t58 * pkin(2) + t59 * qJ(3)) - (pkin(2) * t75 + qJ(3) * t72) * t94) * MDP(14) + (-g(1) * (-t60 * t90 + t61 * t71) - g(2) * (-t58 * t90 + t59 * t71) - (t66 * t85 + t71 * t72) * t94) * MDP(27) + (-g(1) * (t60 * t91 + t61 * t74) - g(2) * (t58 * t91 + t59 * t74) - (-t66 * t86 + t72 * t74) * t94) * MDP(28) + t95 * (g(1) * t61 + g(2) * t59 + g(3) * t89) + (-t70 * MDP(11) + MDP(12) * t68 - t66 * MDP(20) + MDP(21) * t65 - MDP(9)) * t78; t78 * MDP(14); (g(1) * t55 + g(2) * t52 + g(3) * t57) * MDP(21) + (-MDP(27) * t74 + MDP(28) * t71 - MDP(20)) * (g(1) * t54 - g(2) * t80 + g(3) * (-t65 * t89 + t84 * t66)); (-g(1) * t49 + g(2) * t97 - g(3) * (-t57 * t71 - t69 * t85)) * MDP(27) + (g(1) * t50 + g(2) * t96 - g(3) * (-t57 * t74 + t69 * t86)) * MDP(28);];
taug = t1;
