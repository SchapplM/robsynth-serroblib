% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:47
% EndTime: 2022-01-23 08:59:48
% DurationCPUTime: 0.38s
% Computational Cost: add. (122->73), mult. (264->122), div. (0->0), fcn. (286->10), ass. (0->41)
t76 = cos(pkin(9));
t77 = cos(pkin(8));
t79 = sin(qJ(5));
t74 = sin(pkin(8));
t81 = cos(qJ(5));
t95 = t74 * t81;
t58 = t76 * t95 - t79 * t77;
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t73 = sin(pkin(9));
t75 = sin(pkin(7));
t78 = cos(pkin(7));
t92 = t76 * t77;
t56 = t75 * t73 + t78 * t92;
t96 = t74 * t79;
t83 = t56 * t81 + t78 * t96;
t100 = -t82 * t58 + t80 * t83;
t50 = -t56 * t79 + t78 * t95;
t57 = t76 * t96 + t81 * t77;
t99 = -t50 * t82 + t80 * t57;
t98 = t75 * qJ(3) + pkin(1);
t94 = t75 * t80;
t93 = t75 * t82;
t90 = t80 * t74;
t89 = t80 * t77;
t88 = t82 * t74;
t87 = t82 * t77;
t86 = MDP(11) + MDP(15);
t84 = g(1) * t82 + g(2) * t80;
t67 = g(1) * t80 - g(2) * t82;
t71 = t82 * qJ(2);
t70 = t80 * qJ(2);
t65 = pkin(2) * t78 + t98;
t64 = -t74 * pkin(3) + qJ(4) * t77 - qJ(2);
t62 = t78 * t87 + t90;
t61 = t78 * t88 - t89;
t60 = -t78 * t89 + t88;
t59 = t78 * t90 + t87;
t55 = -t78 * t73 + t75 * t92;
t52 = (t77 * pkin(3) + t74 * qJ(4) + pkin(2)) * t78 + t98;
t1 = [(-g(1) * (-t80 * pkin(1) + t71) - g(2) * (t82 * pkin(1) + t70)) * MDP(7) + (-g(1) * t60 - g(2) * t62) * MDP(8) + (-g(1) * (-t65 * t80 + t71) - g(2) * (t65 * t82 + t70)) * MDP(11) + (-g(1) * (t60 * t76 - t73 * t94) - g(2) * (t62 * t76 + t73 * t93)) * MDP(12) + (-g(1) * (-t60 * t73 - t76 * t94) - g(2) * (-t62 * t73 + t76 * t93)) * MDP(13) + (-g(1) * (-t52 * t80 - t64 * t82) - g(2) * (t52 * t82 - t64 * t80)) * MDP(15) + (g(1) * t100 - g(2) * ((t56 * t82 + t76 * t90) * t81 + t61 * t79)) * MDP(21) + (-g(1) * (-t50 * t80 - t82 * t57) + g(2) * t99) * MDP(22) + (MDP(3) - MDP(6)) * t84 + (-MDP(9) + MDP(14)) * (g(1) * t59 - g(2) * t61) + (MDP(2) + t78 * MDP(4) + (-MDP(5) + MDP(10)) * t75) * t67; (-MDP(7) - t86) * t67; t86 * (g(3) * t78 - t75 * t84); (-g(3) * t75 * t74 - g(1) * t61 - g(2) * t59) * MDP(15); (g(1) * t99 - g(2) * (-(t56 * t80 - t76 * t88) * t79 + t59 * t81) - g(3) * (-t79 * t55 + t75 * t95)) * MDP(21) + (-g(1) * (-t80 * t58 - t82 * t83) + g(2) * t100 - g(3) * (-t55 * t81 - t75 * t96)) * MDP(22);];
taug = t1;
