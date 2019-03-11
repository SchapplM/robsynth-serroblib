% Calculate Gravitation load on the joints for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:47
% EndTime: 2019-03-09 03:22:48
% DurationCPUTime: 0.33s
% Computational Cost: add. (150->64), mult. (221->84), div. (0->0), fcn. (185->8), ass. (0->34)
t67 = -qJ(4) - pkin(7);
t70 = sin(qJ(1));
t73 = cos(qJ(1));
t69 = sin(qJ(3));
t87 = t69 * pkin(3);
t97 = t70 * t67 + t73 * t87;
t96 = -g(1) * t70 + g(2) * t73;
t65 = qJ(3) + pkin(9);
t59 = sin(t65);
t60 = cos(t65);
t75 = -g(3) * t59 - t60 * t96;
t68 = sin(qJ(5));
t94 = pkin(5) * t68;
t88 = g(3) * t60;
t86 = t70 * t68;
t71 = cos(qJ(5));
t85 = t70 * t71;
t84 = t73 * t68;
t83 = t73 * t71;
t82 = t73 * pkin(1) + t70 * qJ(2);
t81 = -MDP(15) - MDP(24);
t79 = t70 * t87 + t82;
t62 = t73 * qJ(2);
t78 = -t70 * pkin(1) + t62;
t54 = g(1) * t73 + g(2) * t70;
t58 = t71 * pkin(5) + pkin(4);
t66 = -qJ(6) - pkin(8);
t77 = t59 * t58 + t60 * t66;
t51 = t59 * t84 + t85;
t49 = -t59 * t86 + t83;
t72 = cos(qJ(3));
t52 = t59 * t83 - t86;
t50 = t59 * t85 + t84;
t1 = [(-g(1) * t78 - g(2) * t82) * MDP(6) + (-g(1) * (t78 + t97) - g(2) * (-t73 * t67 + t79)) * MDP(15) + (-g(1) * t52 - g(2) * t50) * MDP(21) + (g(1) * t51 - g(2) * t49) * MDP(22) + (-g(1) * (t62 + t97) - g(2) * t79 + (-g(1) * t77 - g(2) * (-t67 + t94)) * t73 + (-g(1) * (-pkin(1) - t94) - g(2) * t77) * t70) * MDP(24) - (MDP(2) - MDP(4) + MDP(14)) * t96 + (-t69 * MDP(12) - t72 * MDP(13) + t60 * MDP(23) + MDP(3) - MDP(5)) * t54; -(-MDP(6) + t81) * t96; (g(3) * t72 - t69 * t96) * MDP(13) + (t59 * t96 - t88) * MDP(23) + (-g(3) * (-t77 - t87) + t96 * (pkin(3) * t72 + t58 * t60 - t59 * t66)) * MDP(24) + (-MDP(21) * t71 + MDP(22) * t68) * t75 + (pkin(3) * MDP(15) + MDP(12)) * (g(3) * t69 + t96 * t72); t81 * t54; (g(1) * t50 - g(2) * t52 + t71 * t88) * MDP(22) + (pkin(5) * MDP(24) + MDP(21)) * (-g(1) * t49 - g(2) * t51 + t68 * t88); t75 * MDP(24);];
taug  = t1;
