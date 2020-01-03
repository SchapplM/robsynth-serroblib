% Calculate Gravitation load on the joints for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:27
% DurationCPUTime: 0.19s
% Computational Cost: add. (176->53), mult. (231->74), div. (0->0), fcn. (200->8), ass. (0->33)
t72 = cos(qJ(2));
t67 = qJ(2) + qJ(3);
t64 = sin(t67);
t65 = cos(t67);
t81 = t65 * pkin(3) + t64 * qJ(4);
t91 = t72 * pkin(2) + t81;
t90 = MDP(16) - MDP(19);
t89 = MDP(17) - MDP(20);
t88 = pkin(3) * t64;
t86 = g(3) * t65;
t68 = sin(qJ(5));
t70 = sin(qJ(1));
t85 = t70 * t68;
t71 = cos(qJ(5));
t84 = t70 * t71;
t73 = cos(qJ(1));
t83 = t73 * t68;
t82 = t73 * t71;
t80 = qJ(4) * t65;
t69 = sin(qJ(2));
t79 = -pkin(2) * t69 - t88;
t57 = g(1) * t73 + g(2) * t70;
t49 = t57 * t64 - t86;
t77 = t90 * t49 + (-MDP(27) * t68 - MDP(28) * t71 + t89) * (g(3) * t64 + t57 * t65);
t76 = pkin(1) + t91;
t74 = -pkin(7) - pkin(6);
t59 = t73 * t80;
t58 = t70 * t80;
t56 = -t64 * t85 + t82;
t55 = t64 * t84 + t83;
t54 = t64 * t83 + t84;
t53 = t64 * t82 - t85;
t1 = [((g(1) * t74 - g(2) * t76) * t73 + (g(1) * t76 + g(2) * t74) * t70) * MDP(21) + (-g(1) * t56 - g(2) * t54) * MDP(27) + (g(1) * t55 - g(2) * t53) * MDP(28) + (MDP(3) - MDP(18)) * t57 + (-t69 * MDP(10) + t72 * MDP(9) - t89 * t64 + t90 * t65 + MDP(2)) * (g(1) * t70 - g(2) * t73); (-g(3) * t72 + t57 * t69) * MDP(9) + (g(3) * t69 + t57 * t72) * MDP(10) + (-g(1) * (t79 * t73 + t59) - g(2) * (t79 * t70 + t58) - g(3) * t91) * MDP(21) + t77; (-g(1) * (-t73 * t88 + t59) - g(2) * (-t70 * t88 + t58) - g(3) * t81) * MDP(21) + t77; -t49 * MDP(21); (-g(1) * t53 - g(2) * t55 + t71 * t86) * MDP(27) + (g(1) * t54 - g(2) * t56 - t68 * t86) * MDP(28);];
taug = t1;
