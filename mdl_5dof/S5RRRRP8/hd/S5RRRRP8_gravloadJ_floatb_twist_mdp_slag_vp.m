% Calculate Gravitation load on the joints for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:21
% EndTime: 2021-01-16 00:21:23
% DurationCPUTime: 0.32s
% Computational Cost: add. (247->64), mult. (345->93), div. (0->0), fcn. (331->8), ass. (0->38)
t109 = MDP(23) + MDP(25);
t108 = MDP(24) + MDP(26);
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t93 = g(1) * t87 + g(2) * t84;
t110 = MDP(10) - MDP(27);
t83 = sin(qJ(2));
t103 = g(3) * t83;
t86 = cos(qJ(2));
t100 = t84 * t86;
t81 = qJ(3) + qJ(4);
t77 = sin(t81);
t78 = cos(t81);
t98 = t87 * t78;
t62 = t77 * t100 + t98;
t99 = t87 * t77;
t64 = t84 * t78 - t86 * t99;
t56 = -g(1) * t64 + g(2) * t62 + t77 * t103;
t66 = -g(3) * t86 + t93 * t83;
t82 = sin(qJ(3));
t74 = t82 * pkin(3) + pkin(4) * t77;
t101 = pkin(6) + t74;
t97 = t87 * t82;
t85 = cos(qJ(3));
t96 = t87 * t85;
t75 = t85 * pkin(3) + pkin(4) * t78;
t63 = -t78 * t100 + t99;
t65 = t84 * t77 + t86 * t98;
t94 = t108 * (g(1) * t65 - g(2) * t63 + t78 * t103) + t109 * t56;
t73 = pkin(2) + t75;
t80 = -qJ(5) - pkin(8) - pkin(7);
t91 = t86 * t73 - t83 * t80;
t89 = pkin(1) + t91;
t71 = t84 * t82 + t86 * t96;
t70 = t84 * t85 - t86 * t97;
t69 = -t85 * t100 + t97;
t68 = t82 * t100 + t96;
t1 = [t93 * MDP(3) + (-g(1) * t69 - g(2) * t71) * MDP(16) + (-g(1) * t68 - g(2) * t70) * MDP(17) + ((-g(1) * t101 - g(2) * t89) * t87 + (g(1) * t89 - g(2) * t101) * t84) * MDP(28) + t109 * (-g(1) * t63 - g(2) * t65) + t108 * (-g(1) * t62 - g(2) * t64) + (t86 * MDP(9) - t110 * t83 + MDP(2)) * (g(1) * t84 - g(2) * t87); (-g(3) * t91 + t93 * (t73 * t83 + t80 * t86)) * MDP(28) + t110 * (t93 * t86 + t103) + (MDP(16) * t85 - MDP(17) * t82 - t108 * t77 + t109 * t78 + MDP(9)) * t66; (-g(1) * t70 + g(2) * t68 + t82 * t103) * MDP(16) + (g(1) * t71 - g(2) * t69 + t85 * t103) * MDP(17) + (-g(1) * (-t87 * t86 * t74 + t84 * t75) - g(2) * (-t74 * t100 - t87 * t75) + t74 * t103) * MDP(28) + t94; MDP(28) * pkin(4) * t56 + t94; -t66 * MDP(28);];
taug = t1;
