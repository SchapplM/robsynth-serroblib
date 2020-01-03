% Calculate Gravitation load on the joints for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR16_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR16_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:26
% EndTime: 2019-12-31 20:47:28
% DurationCPUTime: 0.39s
% Computational Cost: add. (184->79), mult. (459->137), div. (0->0), fcn. (534->10), ass. (0->37)
t96 = MDP(9) - MDP(12);
t95 = -MDP(10) + MDP(13);
t68 = sin(qJ(2));
t69 = sin(qJ(1));
t72 = cos(qJ(2));
t73 = cos(qJ(1));
t82 = cos(pkin(5));
t80 = t73 * t82;
t60 = t68 * t80 + t69 * t72;
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t59 = t69 * t68 - t72 * t80;
t67 = sin(qJ(4));
t71 = cos(qJ(4));
t65 = sin(pkin(5));
t87 = t65 * t73;
t76 = -t59 * t67 + t71 * t87;
t94 = t60 * t70 + t66 * t76;
t93 = -t60 * t66 + t70 * t76;
t92 = g(3) * t65;
t89 = t65 * t69;
t88 = t65 * t72;
t86 = t66 * t67;
t85 = t66 * t68;
t84 = t67 * t70;
t83 = t68 * t70;
t81 = t69 * t82;
t53 = t59 * t71 + t67 * t87;
t61 = t73 * t68 + t72 * t81;
t50 = -g(1) * t61 - g(2) * t59 + g(3) * t88;
t62 = -t68 * t81 + t73 * t72;
t58 = -t67 * t88 + t71 * t82;
t52 = t61 * t67 + t71 * t89;
t51 = t61 * t71 - t67 * t89;
t49 = t52 * t70 + t62 * t66;
t48 = -t52 * t66 + t62 * t70;
t1 = [(g(1) * t69 - g(2) * t73) * MDP(2) + (-g(1) * (-t69 * pkin(1) - t60 * pkin(2) + pkin(7) * t87 - t59 * qJ(3)) - g(2) * (t73 * pkin(1) + t62 * pkin(2) + pkin(7) * t89 + t61 * qJ(3))) * MDP(14) + (-g(1) * t76 - g(2) * t52) * MDP(20) + (g(1) * t53 - g(2) * t51) * MDP(21) + (-g(1) * t93 - g(2) * t49) * MDP(27) + (g(1) * t94 - g(2) * t48) * MDP(28) + t95 * (g(1) * t59 - g(2) * t61) + t96 * (g(1) * t60 - g(2) * t62) + (-MDP(11) * t65 + MDP(3)) * (g(1) * t73 + g(2) * t69); (-g(1) * (-t61 * pkin(2) + t62 * qJ(3)) - g(2) * (-t59 * pkin(2) + t60 * qJ(3)) - (pkin(2) * t72 + qJ(3) * t68) * t92) * MDP(14) + (-g(1) * (-t61 * t66 + t62 * t84) - g(2) * (-t59 * t66 + t60 * t84) - (t66 * t72 + t67 * t83) * t92) * MDP(27) + (-g(1) * (-t61 * t70 - t62 * t86) - g(2) * (-t59 * t70 - t60 * t86) - (-t67 * t85 + t70 * t72) * t92) * MDP(28) - t96 * t50 + (-MDP(20) * t67 - t71 * MDP(21) - t95) * (g(1) * t62 + g(2) * t60 + t68 * t92); t50 * MDP(14); (g(1) * t52 - g(2) * t76 + g(3) * t58) * MDP(21) + (-MDP(27) * t70 + MDP(28) * t66 - MDP(20)) * (g(1) * t51 + g(2) * t53 + g(3) * (-t67 * t82 - t71 * t88)); (-g(1) * t48 - g(2) * t94 - g(3) * (-t58 * t66 + t65 * t83)) * MDP(27) + (g(1) * t49 - g(2) * t93 - g(3) * (-t58 * t70 - t65 * t85)) * MDP(28);];
taug = t1;
