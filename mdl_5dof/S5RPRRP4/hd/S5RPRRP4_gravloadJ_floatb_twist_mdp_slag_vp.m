% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:52
% EndTime: 2022-01-23 09:32:53
% DurationCPUTime: 0.25s
% Computational Cost: add. (195->61), mult. (264->89), div. (0->0), fcn. (258->8), ass. (0->39)
t107 = MDP(19) + MDP(21);
t106 = MDP(20) + MDP(22);
t81 = qJ(3) + qJ(4);
t74 = cos(t81);
t86 = cos(qJ(3));
t69 = t86 * pkin(3) + pkin(4) * t74;
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t105 = -(-qJ(5) - pkin(7) - pkin(6)) * t82 + (pkin(2) + t69) * t83;
t102 = g(3) * t82;
t87 = cos(qJ(1));
t93 = t87 * t74;
t73 = sin(t81);
t85 = sin(qJ(1));
t99 = t85 * t73;
t59 = t83 * t99 + t93;
t94 = t87 * t73;
t98 = t85 * t74;
t61 = -t83 * t94 + t98;
t55 = -g(1) * t61 + g(2) * t59 + t73 * t102;
t84 = sin(qJ(3));
t97 = t85 * t84;
t96 = t85 * t86;
t68 = t84 * pkin(3) + pkin(4) * t73;
t95 = t87 * t68;
t92 = t87 * t84;
t91 = t87 * t86;
t90 = t87 * pkin(1) + t85 * qJ(2);
t60 = -t83 * t98 + t94;
t62 = t83 * t93 + t99;
t88 = t106 * (g(1) * t62 - g(2) * t60 + t74 * t102) + t107 * t55;
t71 = g(1) * t87 + g(2) * t85;
t70 = g(1) * t85 - g(2) * t87;
t76 = t87 * qJ(2);
t66 = t83 * t91 + t97;
t65 = -t83 * t92 + t96;
t64 = -t83 * t96 + t92;
t63 = t83 * t97 + t91;
t1 = [(-g(1) * (-t85 * pkin(1) + t76) - g(2) * t90) * MDP(6) + (-g(1) * t64 - g(2) * t66) * MDP(12) + (-g(1) * t63 - g(2) * t65) * MDP(13) + (-g(1) * (t76 + t95) - g(2) * (t105 * t87 + t90) + (-g(1) * (-pkin(1) - t105) - g(2) * t68) * t85) * MDP(24) + (MDP(3) - MDP(5)) * t71 + t107 * (-g(1) * t60 - g(2) * t62) + t106 * (-g(1) * t59 - g(2) * t61) + (t82 * MDP(23) + t83 * MDP(4) + MDP(2)) * t70; (-MDP(24) - MDP(6)) * t70; (-g(1) * t65 + g(2) * t63 + t84 * t102) * MDP(12) + (g(1) * t66 - g(2) * t64 + t86 * t102) * MDP(13) + (-g(1) * (t85 * t69 - t83 * t95) - g(2) * (-t85 * t83 * t68 - t87 * t69) + t68 * t102) * MDP(24) + t88; t55 * pkin(4) * MDP(24) + t88; (g(3) * t83 - t71 * t82) * MDP(24);];
taug = t1;
