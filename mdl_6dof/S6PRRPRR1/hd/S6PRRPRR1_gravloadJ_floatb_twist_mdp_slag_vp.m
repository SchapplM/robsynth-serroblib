% Calculate Gravitation load on the joints for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:24
% EndTime: 2021-01-16 03:28:27
% DurationCPUTime: 0.51s
% Computational Cost: add. (374->89), mult. (540->147), div. (0->0), fcn. (624->14), ass. (0->41)
t79 = sin(pkin(6));
t108 = g(3) * t79;
t107 = cos(qJ(2));
t77 = qJ(3) + pkin(12);
t76 = qJ(5) + t77;
t72 = cos(t76);
t82 = sin(qJ(6));
t106 = t72 * t82;
t85 = cos(qJ(6));
t105 = t72 * t85;
t78 = sin(pkin(11));
t104 = t78 * t79;
t84 = sin(qJ(2));
t103 = t78 * t84;
t102 = t79 * t84;
t86 = cos(qJ(3));
t101 = t79 * t86;
t100 = cos(pkin(11));
t99 = t78 * t107;
t98 = t82 * t107;
t97 = t85 * t107;
t96 = t79 * t100;
t95 = t100 * t84;
t80 = cos(pkin(6));
t67 = t80 * t95 + t99;
t71 = sin(t76);
t60 = t67 * t72 - t71 * t96;
t93 = t100 * t107;
t65 = t80 * t103 - t93;
t62 = t71 * t104 - t65 * t72;
t64 = t72 * t102 + t71 * t80;
t94 = (g(1) * t62 + g(2) * t60 + g(3) * t64) * MDP(22) + (-t85 * MDP(28) + t82 * MDP(29) - MDP(21)) * (g(1) * (t72 * t104 + t65 * t71) + g(2) * (-t67 * t71 - t72 * t96) + g(3) * (-t71 * t102 + t72 * t80));
t66 = -t80 * t93 + t103;
t68 = t80 * t99 + t95;
t87 = g(1) * t68 + g(2) * t66 - t107 * t108;
t83 = sin(qJ(3));
t81 = qJ(4) + pkin(8);
t75 = cos(t77);
t74 = sin(t77);
t73 = pkin(3) * t86 + pkin(2);
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(1) * (-t65 * t81 - t68 * t73) - g(2) * (-t66 * t73 + t67 * t81) - (t107 * t73 + t81 * t84) * t108) * MDP(15) + (-g(1) * (-t68 * t105 - t65 * t82) - g(2) * (-t66 * t105 + t67 * t82) - (t72 * t97 + t82 * t84) * t108) * MDP(28) + (-g(1) * (t68 * t106 - t65 * t85) - g(2) * (t66 * t106 + t67 * t85) - (-t72 * t98 + t84 * t85) * t108) * MDP(29) + (MDP(4) - MDP(14)) * (-g(1) * t65 + g(2) * t67 + g(3) * t102) + (t86 * MDP(10) - t83 * MDP(11) + t75 * MDP(12) - t74 * MDP(13) + t72 * MDP(21) - t71 * MDP(22) + MDP(3)) * t87; (-g(1) * (-t83 * t104 + t65 * t86) - g(2) * (-t67 * t86 + t83 * t96) - g(3) * (-t84 * t101 - t80 * t83)) * MDP(11) + (-g(1) * (t75 * t104 + t65 * t74) - g(2) * (-t67 * t74 - t75 * t96) - g(3) * (-t74 * t102 + t75 * t80)) * MDP(12) + (-g(1) * (-t74 * t104 + t65 * t75) - g(2) * (-t67 * t75 + t74 * t96) - g(3) * (-t75 * t102 - t74 * t80)) * MDP(13) + t94 + (pkin(3) * MDP(15) + MDP(10)) * (-g(1) * (t78 * t101 + t65 * t83) - g(2) * (-t67 * t83 - t86 * t96) - g(3) * (-t83 * t102 + t80 * t86)); -t87 * MDP(15); t94; (-g(1) * (-t62 * t82 + t68 * t85) - g(2) * (-t60 * t82 + t66 * t85) - g(3) * (-t64 * t82 - t79 * t97)) * MDP(28) + (-g(1) * (-t62 * t85 - t68 * t82) - g(2) * (-t60 * t85 - t66 * t82) - g(3) * (-t64 * t85 + t79 * t98)) * MDP(29);];
taug = t1;
