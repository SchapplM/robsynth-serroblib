% Calculate Gravitation load on the joints for
% S5RRRRP6
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
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:24
% EndTime: 2021-01-16 00:10:25
% DurationCPUTime: 0.22s
% Computational Cost: add. (236->52), mult. (315->69), div. (0->0), fcn. (284->8), ass. (0->32)
t81 = cos(qJ(2));
t80 = cos(qJ(4));
t70 = t80 * pkin(4) + pkin(3);
t75 = qJ(2) + qJ(3);
t72 = sin(t75);
t73 = cos(t75);
t76 = -qJ(5) - pkin(8);
t90 = t73 * t70 - t72 * t76;
t108 = t81 * pkin(2) + t90;
t107 = MDP(17) - MDP(27);
t106 = MDP(23) + MDP(25);
t105 = MDP(24) + MDP(26);
t79 = sin(qJ(1));
t82 = cos(qJ(1));
t89 = g(1) * t82 + g(2) * t79;
t62 = -g(3) * t73 + t89 * t72;
t99 = g(3) * t72;
t77 = sin(qJ(4));
t96 = t79 * t77;
t95 = t79 * t80;
t94 = t82 * t77;
t93 = t82 * t80;
t91 = pkin(4) * t77 + pkin(6) + pkin(7);
t87 = t70 * t72 + t73 * t76;
t67 = -t73 * t94 + t95;
t65 = t73 * t96 + t93;
t86 = pkin(1) + t108;
t85 = t107 * (t89 * t73 + t99) + (-t105 * t77 + t106 * t80 + MDP(16)) * t62;
t78 = sin(qJ(2));
t68 = t73 * t93 + t96;
t66 = -t73 * t95 + t94;
t1 = [t89 * MDP(3) + ((-g(1) * t91 - g(2) * t86) * t82 + (g(1) * t86 - g(2) * t91) * t79) * MDP(28) + t106 * (-g(1) * t66 - g(2) * t68) + t105 * (-g(1) * t65 - g(2) * t67) + (-t78 * MDP(10) + t73 * MDP(16) + t81 * MDP(9) - t107 * t72 + MDP(2)) * (g(1) * t79 - g(2) * t82); (-g(3) * t81 + t89 * t78) * MDP(9) + (g(3) * t78 + t89 * t81) * MDP(10) + (-g(3) * t108 + t89 * (pkin(2) * t78 + t87)) * MDP(28) + t85; (-g(3) * t90 + t89 * t87) * MDP(28) + t85; t105 * (g(1) * t68 - g(2) * t66 + t80 * t99) + (pkin(4) * MDP(28) + t106) * (-g(1) * t67 + g(2) * t65 + t77 * t99); -t62 * MDP(28);];
taug = t1;
