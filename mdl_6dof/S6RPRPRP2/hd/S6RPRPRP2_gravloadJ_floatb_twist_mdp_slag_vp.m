% Calculate Gravitation load on the joints for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:13
% EndTime: 2019-03-09 03:06:14
% DurationCPUTime: 0.37s
% Computational Cost: add. (305->73), mult. (297->98), div. (0->0), fcn. (275->10), ass. (0->36)
t100 = MDP(19) + MDP(21);
t99 = MDP(20) - MDP(23);
t69 = qJ(3) + pkin(10);
t65 = cos(t69);
t63 = sin(t69);
t95 = t63 * pkin(8);
t101 = t65 * pkin(4) + t95;
t70 = qJ(1) + pkin(9);
t64 = sin(t70);
t66 = cos(t70);
t86 = g(1) * t66 + g(2) * t64;
t96 = g(3) * t63;
t72 = sin(qJ(5));
t93 = t64 * t72;
t75 = cos(qJ(5));
t92 = t64 * t75;
t91 = t66 * t72;
t90 = t66 * t75;
t76 = cos(qJ(3));
t67 = t76 * pkin(3);
t62 = t67 + pkin(2);
t77 = cos(qJ(1));
t89 = t77 * pkin(1) + t66 * t62;
t88 = MDP(13) + MDP(24);
t85 = g(1) * t64 - g(2) * t66;
t71 = -qJ(4) - pkin(7);
t74 = sin(qJ(1));
t83 = -t74 * pkin(1) - t66 * t71;
t82 = pkin(5) * t75 + qJ(6) * t72 + pkin(4);
t54 = t65 * t93 + t90;
t56 = t65 * t91 - t92;
t50 = g(1) * t56 + g(2) * t54 + t72 * t96;
t73 = sin(qJ(3));
t57 = t65 * t90 + t93;
t55 = t65 * t92 - t91;
t1 = [(g(1) * t77 + g(2) * t74) * MDP(3) - t86 * MDP(12) + (-g(1) * (-t64 * t62 + t83) - g(2) * (-t64 * t71 + t89)) * MDP(13) + (-g(1) * (-t55 * pkin(5) - t54 * qJ(6) + t83) - g(2) * (t57 * pkin(5) + t56 * qJ(6) + t101 * t66 + t89) + (-g(1) * (-t62 - t101) + g(2) * t71) * t64) * MDP(24) - t99 * (g(1) * t54 - g(2) * t56) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t74 - g(2) * t77) + t100 * (g(1) * t55 - g(2) * t57) + (t76 * MDP(10) - t73 * MDP(11) + t63 * MDP(22)) * t85; (-MDP(4) - t88) * g(3); (g(3) * t73 + t86 * t76) * MDP(11) + (-t86 * t65 - t96) * MDP(22) + (-g(3) * (t82 * t65 + t67 + t95) + t86 * (pkin(3) * t73 - pkin(8) * t65 + t82 * t63)) * MDP(24) + (pkin(3) * MDP(13) + MDP(10)) * (-g(3) * t76 + t86 * t73) + (t100 * t75 - t99 * t72) * (-g(3) * t65 + t86 * t63); -t88 * t85; (-g(1) * (-t56 * pkin(5) + t57 * qJ(6)) - g(2) * (-t54 * pkin(5) + t55 * qJ(6)) - (-pkin(5) * t72 + qJ(6) * t75) * t96) * MDP(24) + t99 * (g(1) * t57 + g(2) * t55 + t75 * t96) + t100 * t50; -t50 * MDP(24);];
taug  = t1;
