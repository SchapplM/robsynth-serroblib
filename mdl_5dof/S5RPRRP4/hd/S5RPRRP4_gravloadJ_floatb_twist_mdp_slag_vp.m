% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:37
% EndTime: 2021-01-15 12:56:39
% DurationCPUTime: 0.26s
% Computational Cost: add. (195->61), mult. (264->89), div. (0->0), fcn. (258->8), ass. (0->39)
t105 = MDP(19) + MDP(21);
t104 = MDP(20) + MDP(22);
t79 = qJ(3) + qJ(4);
t73 = cos(t79);
t84 = cos(qJ(3));
t68 = t84 * pkin(3) + pkin(4) * t73;
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t103 = -(-qJ(5) - pkin(7) - pkin(6)) * t80 + (pkin(2) + t68) * t81;
t102 = g(1) * t80;
t85 = cos(qJ(1));
t91 = t85 * t73;
t72 = sin(t79);
t83 = sin(qJ(1));
t96 = t83 * t72;
t58 = -t81 * t96 - t91;
t92 = t85 * t72;
t95 = t83 * t73;
t60 = t81 * t92 - t95;
t54 = -g(2) * t58 - g(3) * t60 + t72 * t102;
t82 = sin(qJ(3));
t67 = t82 * pkin(3) + pkin(4) * t72;
t97 = t83 * t67;
t94 = t83 * t82;
t93 = t83 * t84;
t90 = t85 * t82;
t89 = t85 * t84;
t88 = t85 * pkin(1) + t83 * qJ(2);
t59 = t81 * t95 - t92;
t61 = t81 * t91 + t96;
t86 = t104 * (g(2) * t59 - g(3) * t61 + t73 * t102) + t105 * t54;
t70 = g(2) * t85 + g(3) * t83;
t69 = g(2) * t83 - g(3) * t85;
t75 = t83 * pkin(1);
t65 = t81 * t89 + t94;
t64 = t81 * t90 - t93;
t63 = t81 * t93 - t90;
t62 = -t81 * t94 - t89;
t1 = [(-g(2) * t88 - g(3) * (-t85 * qJ(2) + t75)) * MDP(6) + (-g(2) * t65 - g(3) * t63) * MDP(12) + (g(2) * t64 - g(3) * t62) * MDP(13) + (-g(2) * (t88 + t97) - g(3) * (t103 * t83 + t75) + (-g(2) * t103 - g(3) * (-qJ(2) - t67)) * t85) * MDP(24) + (MDP(3) - MDP(5)) * t69 + t105 * (-g(2) * t61 - g(3) * t59) + t104 * (g(2) * t60 - g(3) * t58) + (-t80 * MDP(23) - t81 * MDP(4) - MDP(2)) * t70; (MDP(24) + MDP(6)) * t70; (-g(2) * t62 - g(3) * t64 + t82 * t102) * MDP(12) + (g(2) * t63 - g(3) * t65 + t84 * t102) * MDP(13) + (t67 * t102 - g(2) * (-t85 * t68 - t81 * t97) - g(3) * (t85 * t81 * t67 - t83 * t68)) * MDP(24) + t86; t54 * pkin(4) * MDP(24) + t86; (g(1) * t81 - t69 * t80) * MDP(24);];
taug = t1;
