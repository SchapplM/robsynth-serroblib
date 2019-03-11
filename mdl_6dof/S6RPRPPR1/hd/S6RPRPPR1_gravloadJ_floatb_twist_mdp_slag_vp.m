% Calculate Gravitation load on the joints for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:35
% EndTime: 2019-03-09 02:39:36
% DurationCPUTime: 0.38s
% Computational Cost: add. (241->69), mult. (214->98), div. (0->0), fcn. (184->12), ass. (0->38)
t64 = qJ(3) + pkin(10);
t56 = sin(t64);
t59 = cos(t64);
t97 = t59 * pkin(4) + t56 * qJ(5);
t65 = qJ(1) + pkin(9);
t57 = sin(t65);
t60 = cos(t65);
t81 = g(1) * t60 + g(2) * t57;
t75 = -g(3) * t59 + t81 * t56;
t94 = g(3) * t56;
t91 = t57 * t59;
t66 = sin(pkin(11));
t90 = t57 * t66;
t67 = cos(pkin(11));
t89 = t57 * t67;
t63 = pkin(11) + qJ(6);
t55 = sin(t63);
t88 = t60 * t55;
t58 = cos(t63);
t87 = t60 * t58;
t86 = t60 * t66;
t85 = t60 * t67;
t71 = cos(qJ(3));
t61 = t71 * pkin(3);
t54 = t61 + pkin(2);
t72 = cos(qJ(1));
t84 = t72 * pkin(1) + t60 * t54;
t82 = MDP(13) + MDP(17);
t80 = g(1) * t57 - g(2) * t60;
t68 = -qJ(4) - pkin(7);
t70 = sin(qJ(1));
t78 = -t70 * pkin(1) - t60 * t68;
t69 = sin(qJ(3));
t51 = t57 * t55 + t59 * t87;
t50 = t57 * t58 - t59 * t88;
t49 = -t58 * t91 + t88;
t48 = t55 * t91 + t87;
t1 = [(g(1) * t72 + g(2) * t70) * MDP(3) - t81 * MDP(12) + (-g(1) * (-t57 * t54 + t78) - g(2) * (-t57 * t68 + t84)) * MDP(13) + (-g(1) * (-t59 * t89 + t86) - g(2) * (t59 * t85 + t90)) * MDP(14) + (-g(1) * (t59 * t90 + t85) - g(2) * (-t59 * t86 + t89)) * MDP(15) + (-g(1) * t78 - g(2) * (t60 * t97 + t84) + (-g(1) * (-t54 - t97) + g(2) * t68) * t57) * MDP(17) + (-g(1) * t49 - g(2) * t51) * MDP(23) + (-g(1) * t48 - g(2) * t50) * MDP(24) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t70 - g(2) * t72) + (t71 * MDP(10) - t69 * MDP(11) + t56 * MDP(16)) * t80; (-MDP(4) - t82) * g(3); (g(3) * t69 + t81 * t71) * MDP(11) + (-t81 * t59 - t94) * MDP(16) + (-g(3) * (t61 + t97) + t81 * (pkin(3) * t69 + pkin(4) * t56 - qJ(5) * t59)) * MDP(17) + (MDP(13) * pkin(3) + MDP(10)) * (-g(3) * t71 + t81 * t69) + (MDP(14) * t67 - MDP(15) * t66 + MDP(23) * t58 - MDP(24) * t55) * t75; -t82 * t80; -t75 * MDP(17); (-g(1) * t50 + g(2) * t48 + t55 * t94) * MDP(23) + (g(1) * t51 - g(2) * t49 + t58 * t94) * MDP(24);];
taug  = t1;
