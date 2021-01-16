% Calculate Gravitation load on the joints for
% S5RRRPR9
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
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:57
% EndTime: 2021-01-15 23:23:59
% DurationCPUTime: 0.24s
% Computational Cost: add. (222->72), mult. (287->105), div. (0->0), fcn. (275->10), ass. (0->42)
t106 = MDP(10) - MDP(20);
t83 = sin(qJ(2));
t86 = cos(qJ(2));
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t90 = g(1) * t87 + g(2) * t84;
t63 = -g(3) * t86 + t83 * t90;
t102 = g(3) * t83;
t100 = t84 * t86;
t80 = qJ(3) + pkin(9);
t79 = qJ(5) + t80;
t73 = sin(t79);
t99 = t87 * t73;
t74 = cos(t79);
t98 = t87 * t74;
t77 = sin(t80);
t97 = t87 * t77;
t78 = cos(t80);
t96 = t87 * t78;
t82 = sin(qJ(3));
t95 = t87 * t82;
t85 = cos(qJ(3));
t94 = t87 * t85;
t55 = t73 * t100 + t98;
t56 = -t74 * t100 + t99;
t57 = t84 * t74 - t86 * t99;
t58 = t84 * t73 + t86 * t98;
t93 = (-g(1) * t57 + g(2) * t55 + t73 * t102) * MDP(27) + (g(1) * t58 - g(2) * t56 + t74 * t102) * MDP(28);
t76 = t85 * pkin(3) + pkin(2);
t81 = qJ(4) + pkin(7);
t91 = t76 * t86 + t81 * t83;
t69 = t84 * t85 - t86 * t95;
t67 = t82 * t100 + t94;
t75 = t82 * pkin(3) + pkin(6);
t70 = t84 * t82 + t86 * t94;
t68 = -t85 * t100 + t95;
t65 = pkin(1) + t91;
t62 = t84 * t77 + t86 * t96;
t61 = t84 * t78 - t86 * t97;
t60 = -t78 * t100 + t97;
t59 = t77 * t100 + t96;
t1 = [t90 * MDP(3) + (-g(1) * t68 - g(2) * t70) * MDP(16) + (-g(1) * t67 - g(2) * t69) * MDP(17) + (-g(1) * t60 - g(2) * t62) * MDP(18) + (-g(1) * t59 - g(2) * t61) * MDP(19) + (-g(1) * (-t65 * t84 + t75 * t87) - g(2) * (t65 * t87 + t75 * t84)) * MDP(21) + (-g(1) * t56 - g(2) * t58) * MDP(27) + (-g(1) * t55 - g(2) * t57) * MDP(28) + (MDP(9) * t86 - t106 * t83 + MDP(2)) * (g(1) * t84 - g(2) * t87); (-g(3) * t91 - t90 * (-t83 * t76 + t81 * t86)) * MDP(21) + t106 * (t86 * t90 + t102) + (t85 * MDP(16) - t82 * MDP(17) + t78 * MDP(18) - t77 * MDP(19) + t74 * MDP(27) - t73 * MDP(28) + MDP(9)) * t63; (g(1) * t70 - g(2) * t68 + t85 * t102) * MDP(17) + (-g(1) * t61 + g(2) * t59 + t77 * t102) * MDP(18) + (g(1) * t62 - g(2) * t60 + t78 * t102) * MDP(19) + t93 + (pkin(3) * MDP(21) + MDP(16)) * (-g(1) * t69 + g(2) * t67 + t82 * t102); -t63 * MDP(21); t93;];
taug = t1;
