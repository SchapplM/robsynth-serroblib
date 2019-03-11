% Calculate Gravitation load on the joints for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:53
% EndTime: 2019-03-09 03:52:54
% DurationCPUTime: 0.32s
% Computational Cost: add. (288->73), mult. (276->104), div. (0->0), fcn. (258->12), ass. (0->44)
t85 = sin(qJ(1));
t86 = cos(qJ(1));
t68 = g(1) * t86 + g(2) * t85;
t108 = MDP(14) - MDP(17);
t79 = pkin(10) + qJ(3);
t73 = sin(t79);
t75 = cos(t79);
t56 = -g(3) * t75 + t68 * t73;
t105 = g(3) * t73;
t78 = pkin(11) + qJ(5);
t76 = qJ(6) + t78;
t69 = sin(t76);
t103 = t85 * t69;
t70 = cos(t76);
t102 = t85 * t70;
t72 = sin(t78);
t101 = t85 * t72;
t74 = cos(t78);
t100 = t85 * t74;
t80 = sin(pkin(11));
t99 = t85 * t80;
t82 = cos(pkin(11));
t98 = t85 * t82;
t97 = t86 * t69;
t96 = t86 * t70;
t95 = t86 * t72;
t94 = t86 * t74;
t93 = t86 * t80;
t92 = t86 * t82;
t58 = t103 * t75 + t96;
t59 = -t102 * t75 + t97;
t60 = -t75 * t97 + t102;
t61 = t75 * t96 + t103;
t91 = (-g(1) * t60 + g(2) * t58 + t105 * t69) * MDP(31) + (g(1) * t61 - g(2) * t59 + t105 * t70) * MDP(32);
t67 = g(1) * t85 - g(2) * t86;
t90 = t75 * pkin(3) + t73 * qJ(4);
t83 = cos(pkin(10));
t88 = t83 * pkin(2) + pkin(1) + t90;
t84 = -pkin(7) - qJ(2);
t65 = t75 * t94 + t101;
t64 = -t75 * t95 + t100;
t63 = -t100 * t75 + t95;
t62 = t101 * t75 + t94;
t1 = [(-g(1) * (-t85 * pkin(1) + t86 * qJ(2)) - g(2) * (t86 * pkin(1) + t85 * qJ(2))) * MDP(7) + (-g(1) * (-t75 * t98 + t93) - g(2) * (t75 * t92 + t99)) * MDP(15) + (-g(1) * (t75 * t99 + t92) - g(2) * (-t75 * t93 + t98)) * MDP(16) + ((g(1) * t84 - g(2) * t88) * t86 + (g(1) * t88 + g(2) * t84) * t85) * MDP(18) + (-g(1) * t63 - g(2) * t65) * MDP(24) + (-g(1) * t62 - g(2) * t64) * MDP(25) + (-g(1) * t59 - g(2) * t61) * MDP(31) + (-g(1) * t58 - g(2) * t60) * MDP(32) + (MDP(3) - MDP(6)) * t68 + (t75 * MDP(13) + MDP(4) * t83 - MDP(5) * sin(pkin(10)) - t108 * t73 + MDP(2)) * t67; (-MDP(18) - MDP(7)) * t67; (-g(3) * t90 + t68 * (pkin(3) * t73 - qJ(4) * t75)) * MDP(18) + t108 * (t68 * t75 + t105) + (MDP(15) * t82 - MDP(16) * t80 + MDP(24) * t74 - MDP(25) * t72 + MDP(31) * t70 - MDP(32) * t69 + MDP(13)) * t56; -t56 * MDP(18); (-g(1) * t64 + g(2) * t62 + t105 * t72) * MDP(24) + (g(1) * t65 - g(2) * t63 + t105 * t74) * MDP(25) + t91; t91;];
taug  = t1;
