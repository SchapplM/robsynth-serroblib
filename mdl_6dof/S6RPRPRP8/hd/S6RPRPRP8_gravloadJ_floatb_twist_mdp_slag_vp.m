% Calculate Gravitation load on the joints for
% S6RPRPRP8
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:05
% EndTime: 2019-03-09 03:26:06
% DurationCPUTime: 0.39s
% Computational Cost: add. (212->72), mult. (315->94), div. (0->0), fcn. (291->8), ass. (0->34)
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t105 = -g(1) * t76 + g(2) * t79;
t104 = MDP(21) + MDP(23);
t103 = MDP(22) - MDP(25);
t72 = qJ(3) + pkin(9);
t67 = cos(t72);
t97 = g(3) * t67;
t96 = t67 * pkin(8);
t75 = sin(qJ(3));
t95 = t75 * pkin(3);
t74 = sin(qJ(5));
t94 = t74 * t79;
t93 = t76 * t74;
t77 = cos(qJ(5));
t92 = t76 * t77;
t91 = t79 * t77;
t90 = t79 * pkin(1) + t76 * qJ(2);
t89 = -MDP(15) - MDP(26);
t88 = -t76 * pkin(1) + t79 * qJ(2);
t66 = sin(t72);
t86 = pkin(4) * t66 - t96;
t61 = g(1) * t79 + g(2) * t76;
t73 = -qJ(4) - pkin(7);
t84 = t76 * t73 + t79 * t95 + t88;
t83 = -t73 * t79 + t76 * t95 + t90;
t82 = pkin(5) * t77 + qJ(6) * t74 + pkin(4);
t56 = t66 * t93 - t91;
t58 = t66 * t94 + t92;
t52 = g(1) * t56 - g(2) * t58 + t74 * t97;
t78 = cos(qJ(3));
t59 = t66 * t91 - t93;
t57 = t66 * t92 + t94;
t1 = [(-g(1) * t88 - g(2) * t90) * MDP(6) + (-g(1) * t84 - g(2) * t83) * MDP(15) + (-g(1) * (t59 * pkin(5) + t58 * qJ(6) + t86 * t79 + t84) - g(2) * (t57 * pkin(5) + t56 * qJ(6) + t86 * t76 + t83)) * MDP(26) + t104 * (-g(1) * t59 - g(2) * t57) - (MDP(2) - MDP(4) + MDP(14)) * t105 + t103 * (g(1) * t58 + g(2) * t56) + (-t75 * MDP(12) - t78 * MDP(13) + t67 * MDP(24) + MDP(3) - MDP(5)) * t61; -(-MDP(6) + t89) * t105; (g(3) * t78 - t105 * t75) * MDP(13) + (t105 * t66 - t97) * MDP(24) + (-g(3) * (-t82 * t66 - t95 + t96) + t105 * (pkin(3) * t78 + pkin(8) * t66 + t82 * t67)) * MDP(26) + (pkin(3) * MDP(15) + MDP(12)) * (g(3) * t75 + t105 * t78) + (t103 * t74 - t104 * t77) * (-g(3) * t66 - t105 * t67); t89 * t61; (-g(1) * (-pkin(5) * t56 + qJ(6) * t57) - g(2) * (pkin(5) * t58 - qJ(6) * t59) - (-pkin(5) * t74 + qJ(6) * t77) * t97) * MDP(26) + t103 * (g(1) * t57 - g(2) * t59 + t77 * t97) + t104 * t52; -t52 * MDP(26);];
taug  = t1;
