% Calculate Gravitation load on the joints for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:45
% EndTime: 2019-03-08 20:07:47
% DurationCPUTime: 0.44s
% Computational Cost: add. (310->81), mult. (540->130), div. (0->0), fcn. (616->12), ass. (0->43)
t107 = MDP(15) - MDP(23);
t81 = sin(qJ(2));
t83 = cos(qJ(2));
t97 = cos(pkin(10));
t98 = cos(pkin(6));
t91 = t98 * t97;
t96 = sin(pkin(10));
t61 = t81 * t96 - t83 * t91;
t106 = g(2) * t61;
t62 = t81 * t91 + t83 * t96;
t105 = g(2) * t62;
t76 = sin(pkin(6));
t104 = g(3) * t76;
t74 = pkin(11) + qJ(4);
t73 = cos(t74);
t80 = sin(qJ(5));
t103 = t73 * t80;
t82 = cos(qJ(5));
t102 = t73 * t82;
t101 = t76 * t81;
t100 = t80 * t83;
t99 = t82 * t83;
t95 = MDP(24) + MDP(8);
t94 = pkin(5) * t80 + pkin(8) + qJ(3);
t93 = t76 * t97;
t92 = t76 * t96;
t90 = t98 * t96;
t71 = pkin(5) * t82 + pkin(4);
t72 = sin(t74);
t77 = cos(pkin(11));
t78 = -qJ(6) - pkin(9);
t89 = pkin(3) * t77 + t71 * t73 - t72 * t78 + pkin(2);
t55 = t62 * t72 + t73 * t93;
t64 = -t81 * t90 + t83 * t97;
t57 = t64 * t72 - t73 * t92;
t59 = t101 * t72 - t73 * t98;
t88 = g(1) * t57 + g(2) * t55 + g(3) * t59;
t63 = t81 * t97 + t83 * t90;
t86 = -g(1) * t63 + t104 * t83 - t106;
t60 = t101 * t73 + t72 * t98;
t58 = t64 * t73 + t72 * t92;
t56 = t62 * t73 - t72 * t93;
t1 = [(-MDP(1) - t95) * g(3); (-g(1) * (-pkin(2) * t63 + qJ(3) * t64) - g(2) * (-pkin(2) * t61 + qJ(3) * t62) - (pkin(2) * t83 + qJ(3) * t81) * t104) * MDP(8) + (-g(1) * (-t102 * t63 + t64 * t80) - g(2) * (-t102 * t61 + t62 * t80) - (t73 * t99 + t80 * t81) * t104) * MDP(21) + (-g(1) * (t103 * t63 + t64 * t82) - g(2) * (t103 * t61 + t62 * t82) - (-t100 * t73 + t81 * t82) * t104) * MDP(22) + (-g(1) * (-t63 * t89 + t64 * t94) - t94 * t105 + t89 * t106 - (t81 * t94 + t83 * t89) * t104) * MDP(24) + (MDP(4) - MDP(7)) * (g(1) * t64 + g(3) * t101 + t105) + (-t73 * MDP(14) - t77 * MDP(5) + MDP(6) * sin(pkin(11)) + t107 * t72 - MDP(3)) * t86; t95 * t86; (-g(1) * (-t57 * t71 - t58 * t78) - g(2) * (-t55 * t71 - t56 * t78) - g(3) * (-t59 * t71 - t60 * t78)) * MDP(24) + t107 * (g(1) * t58 + g(2) * t56 + g(3) * t60) + (MDP(21) * t82 - MDP(22) * t80 + MDP(14)) * t88; (-g(1) * (-t58 * t82 - t63 * t80) - g(2) * (-t56 * t82 - t61 * t80) - g(3) * (t100 * t76 - t60 * t82)) * MDP(22) + (MDP(24) * pkin(5) + MDP(21)) * (-g(1) * (-t58 * t80 + t63 * t82) - g(2) * (-t56 * t80 + t61 * t82) - g(3) * (-t60 * t80 - t76 * t99)); -t88 * MDP(24);];
taug  = t1;
