% Calculate Gravitation load on the joints for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:54
% EndTime: 2021-01-15 19:25:56
% DurationCPUTime: 0.23s
% Computational Cost: add. (110->48), mult. (226->64), div. (0->0), fcn. (203->6), ass. (0->26)
t67 = cos(qJ(4));
t61 = t67 * pkin(4) + pkin(3);
t63 = -qJ(5) - pkin(7);
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t91 = -t61 * t65 - t63 * t68;
t90 = qJ(2) - t91;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t89 = -g(1) * t66 + g(2) * t69;
t75 = MDP(19) + MDP(21);
t74 = MDP(20) + MDP(22);
t51 = -g(3) * t65 - t68 * t89;
t82 = g(3) * t68;
t64 = sin(qJ(4));
t80 = t66 * t64;
t79 = t66 * t67;
t78 = t69 * t64;
t77 = t69 * t67;
t76 = MDP(13) - MDP(23);
t54 = t65 * t78 + t79;
t52 = -t65 * t80 + t77;
t59 = t64 * pkin(4) + pkin(1) + pkin(6);
t55 = t65 * t77 - t80;
t53 = t65 * t79 + t78;
t1 = [(-g(1) * (-t66 * pkin(1) + qJ(2) * t69) - g(2) * (pkin(1) * t69 + t66 * qJ(2))) * MDP(6) + (-g(1) * (-t59 * t66 + t90 * t69) - g(2) * (t59 * t69 + t90 * t66)) * MDP(24) - (MDP(2) - MDP(4)) * t89 + t75 * (-g(1) * t55 - g(2) * t53) + t74 * (g(1) * t54 - g(2) * t52) + (-t65 * MDP(12) - t76 * t68 + MDP(3) - MDP(5)) * (g(1) * t69 + g(2) * t66); -(-MDP(24) - MDP(6)) * t89; (-g(3) * t91 + t89 * (t61 * t68 - t63 * t65)) * MDP(24) + t76 * (-t65 * t89 + t82) + (t74 * t64 - t75 * t67 - MDP(12)) * t51; t74 * (g(1) * t53 - g(2) * t55 + t67 * t82) + (pkin(4) * MDP(24) + t75) * (-g(1) * t52 - g(2) * t54 + t64 * t82); t51 * MDP(24);];
taug = t1;
