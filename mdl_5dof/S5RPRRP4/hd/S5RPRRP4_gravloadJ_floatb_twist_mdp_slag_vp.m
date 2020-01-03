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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:24
% EndTime: 2020-01-03 11:50:26
% DurationCPUTime: 0.26s
% Computational Cost: add. (141->60), mult. (207->89), div. (0->0), fcn. (194->8), ass. (0->37)
t73 = qJ(3) + qJ(4);
t67 = cos(t73);
t78 = cos(qJ(3));
t62 = t78 * pkin(3) + pkin(4) * t67;
t74 = sin(pkin(8));
t75 = cos(pkin(8));
t98 = -(-qJ(5) - pkin(7) - pkin(6)) * t74 + (pkin(2) + t62) * t75;
t79 = cos(qJ(1));
t85 = t79 * t67;
t66 = sin(t73);
t77 = sin(qJ(1));
t90 = t77 * t66;
t51 = -t75 * t90 - t85;
t86 = t79 * t66;
t89 = t77 * t67;
t53 = t75 * t86 - t89;
t96 = g(1) * t74;
t97 = -g(2) * t51 - g(3) * t53 + t66 * t96;
t76 = sin(qJ(3));
t61 = t76 * pkin(3) + pkin(4) * t66;
t91 = t77 * t61;
t88 = t77 * t76;
t87 = t77 * t78;
t84 = t79 * t76;
t83 = t79 * t78;
t52 = t75 * t89 - t86;
t54 = t75 * t85 + t90;
t82 = t97 * MDP(20) + (g(2) * t52 - g(3) * t54 + t67 * t96) * MDP(21);
t81 = t79 * pkin(1) + t77 * qJ(2);
t64 = g(2) * t79 + g(3) * t77;
t63 = g(2) * t77 - g(3) * t79;
t69 = t77 * pkin(1);
t58 = t75 * t83 + t88;
t57 = t75 * t84 - t87;
t56 = t75 * t87 - t84;
t55 = -t75 * t88 - t83;
t1 = [(-g(2) * t81 - g(3) * (-t79 * qJ(2) + t69)) * MDP(7) + (-g(2) * t58 - g(3) * t56) * MDP(13) + (g(2) * t57 - g(3) * t55) * MDP(14) + (-g(2) * t54 - g(3) * t52) * MDP(20) + (g(2) * t53 - g(3) * t51) * MDP(21) + (-g(2) * (t81 + t91) - g(3) * (t98 * t77 + t69) + (-g(2) * t98 - g(3) * (-qJ(2) - t61)) * t79) * MDP(23) + (MDP(3) - MDP(6)) * t63 + (-t75 * MDP(4) - MDP(2) + (MDP(5) - MDP(22)) * t74) * t64; (MDP(23) + MDP(7)) * t64; (-g(2) * t55 - g(3) * t57 + t76 * t96) * MDP(13) + (g(2) * t56 - g(3) * t58 + t78 * t96) * MDP(14) + (t61 * t96 - g(2) * (-t79 * t62 - t75 * t91) - g(3) * (t79 * t75 * t61 - t77 * t62)) * MDP(23) + t82; t97 * MDP(23) * pkin(4) + t82; (g(1) * t75 - t63 * t74) * MDP(23);];
taug = t1;
