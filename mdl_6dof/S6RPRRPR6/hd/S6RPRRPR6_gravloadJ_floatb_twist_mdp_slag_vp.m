% Calculate Gravitation load on the joints for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:43
% EndTime: 2019-03-09 05:17:44
% DurationCPUTime: 0.28s
% Computational Cost: add. (260->66), mult. (265->91), div. (0->0), fcn. (243->10), ass. (0->39)
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t68 = g(1) * t85 + g(2) * t83;
t108 = MDP(14) - MDP(22);
t77 = pkin(10) + qJ(3);
t73 = sin(t77);
t74 = cos(t77);
t56 = -g(3) * t74 + t68 * t73;
t102 = g(3) * t73;
t75 = qJ(4) + pkin(11) + qJ(6);
t69 = sin(t75);
t100 = t83 * t69;
t70 = cos(t75);
t99 = t83 * t70;
t82 = sin(qJ(4));
t98 = t83 * t82;
t84 = cos(qJ(4));
t97 = t83 * t84;
t96 = t85 * t69;
t95 = t85 * t70;
t94 = t85 * t82;
t93 = t85 * t84;
t58 = t74 * t100 + t95;
t59 = -t74 * t99 + t96;
t60 = -t74 * t96 + t99;
t61 = t74 * t95 + t100;
t92 = (-g(1) * t60 + g(2) * t58 + t69 * t102) * MDP(29) + (g(1) * t61 - g(2) * t59 + t70 * t102) * MDP(30);
t90 = pkin(4) * t82 + pkin(7) + qJ(2);
t67 = g(1) * t83 - g(2) * t85;
t72 = t84 * pkin(4) + pkin(3);
t80 = -qJ(5) - pkin(8);
t89 = t74 * t72 - t73 * t80;
t65 = -t74 * t94 + t97;
t63 = t74 * t98 + t93;
t79 = cos(pkin(10));
t87 = t79 * pkin(2) + pkin(1) + t89;
t66 = t74 * t93 + t98;
t64 = -t74 * t97 + t94;
t1 = [(-g(1) * (-t83 * pkin(1) + t85 * qJ(2)) - g(2) * (t85 * pkin(1) + t83 * qJ(2))) * MDP(7) + (-g(1) * t64 - g(2) * t66) * MDP(20) + (-g(1) * t63 - g(2) * t65) * MDP(21) + ((-g(1) * t90 - g(2) * t87) * t85 + (g(1) * t87 - g(2) * t90) * t83) * MDP(23) + (-g(1) * t59 - g(2) * t61) * MDP(29) + (-g(1) * t58 - g(2) * t60) * MDP(30) + (MDP(3) - MDP(6)) * t68 + (t74 * MDP(13) + MDP(4) * t79 - MDP(5) * sin(pkin(10)) - t108 * t73 + MDP(2)) * t67; (-MDP(23) - MDP(7)) * t67; (-g(3) * t89 + t68 * (t72 * t73 + t74 * t80)) * MDP(23) + t108 * (t68 * t74 + t102) + (MDP(20) * t84 - MDP(21) * t82 + MDP(29) * t70 - MDP(30) * t69 + MDP(13)) * t56; (g(1) * t66 - g(2) * t64 + t84 * t102) * MDP(21) + t92 + (pkin(4) * MDP(23) + MDP(20)) * (-g(1) * t65 + g(2) * t63 + t82 * t102); -t56 * MDP(23); t92;];
taug  = t1;
