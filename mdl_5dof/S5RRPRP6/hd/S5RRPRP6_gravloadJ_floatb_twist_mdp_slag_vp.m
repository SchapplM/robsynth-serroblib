% Calculate Gravitation load on the joints for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:59
% EndTime: 2021-01-15 20:30:01
% DurationCPUTime: 0.32s
% Computational Cost: add. (193->60), mult. (271->76), div. (0->0), fcn. (241->8), ass. (0->33)
t102 = MDP(12) - MDP(24);
t100 = MDP(20) + MDP(22);
t99 = MDP(21) + MDP(23);
t79 = cos(qJ(4));
t67 = t79 * pkin(4) + pkin(3);
t73 = qJ(2) + pkin(8);
t69 = sin(t73);
t70 = cos(t73);
t74 = -qJ(5) - pkin(7);
t101 = t70 * t67 - t69 * t74;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t64 = g(1) * t81 + g(2) * t78;
t56 = -g(3) * t70 + t64 * t69;
t93 = g(3) * t69;
t76 = sin(qJ(4));
t89 = t78 * t76;
t88 = t78 * t79;
t87 = t81 * t76;
t86 = t81 * t79;
t63 = g(1) * t78 - g(2) * t81;
t61 = -t70 * t87 + t88;
t59 = t70 * t89 + t86;
t80 = cos(qJ(2));
t77 = sin(qJ(2));
t75 = -qJ(3) - pkin(6);
t71 = t80 * pkin(2);
t68 = t71 + pkin(1);
t66 = t75 * t81;
t65 = t81 * t68;
t62 = t70 * t86 + t89;
t60 = -t70 * t88 + t87;
t1 = [(-g(1) * (-t78 * t68 - t66) - g(2) * (-t78 * t75 + t65)) * MDP(14) + (-g(1) * (pkin(4) * t87 - t66) - g(2) * (t101 * t81 + t65) + (-g(1) * (-t68 - t101) - g(2) * (pkin(4) * t76 - t75)) * t78) * MDP(25) + (MDP(3) - MDP(13)) * t64 + t100 * (-g(1) * t60 - g(2) * t62) + t99 * (-g(1) * t59 - g(2) * t61) + (-t77 * MDP(10) + t70 * MDP(11) + t80 * MDP(9) - t102 * t69 + MDP(2)) * t63; (g(3) * t77 + t64 * t80) * MDP(10) + (-g(3) * (t71 + t101) + t64 * (pkin(2) * t77 + t67 * t69 + t70 * t74)) * MDP(25) + (pkin(2) * MDP(14) + MDP(9)) * (-g(3) * t80 + t64 * t77) + t102 * (t64 * t70 + t93) + (t100 * t79 - t99 * t76 + MDP(11)) * t56; (-MDP(14) - MDP(25)) * t63; t99 * (g(1) * t62 - g(2) * t60 + t79 * t93) + (pkin(4) * MDP(25) + t100) * (-g(1) * t61 + g(2) * t59 + t76 * t93); -t56 * MDP(25);];
taug = t1;
