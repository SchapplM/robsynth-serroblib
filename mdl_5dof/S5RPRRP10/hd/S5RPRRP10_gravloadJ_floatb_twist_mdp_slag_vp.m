% Calculate Gravitation load on the joints for
% S5RPRRP10
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
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:50
% EndTime: 2021-01-15 19:14:52
% DurationCPUTime: 0.27s
% Computational Cost: add. (175->50), mult. (238->66), div. (0->0), fcn. (215->7), ass. (0->29)
t69 = sin(qJ(1));
t71 = cos(qJ(1));
t58 = g(1) * t71 + g(2) * t69;
t90 = MDP(13) - MDP(23);
t89 = MDP(19) + MDP(21);
t88 = MDP(20) + MDP(22);
t64 = pkin(8) + qJ(3);
t61 = sin(t64);
t62 = cos(t64);
t50 = -g(3) * t62 + t58 * t61;
t82 = g(3) * t61;
t68 = sin(qJ(4));
t80 = t69 * t68;
t70 = cos(qJ(4));
t79 = t69 * t70;
t78 = t71 * t68;
t77 = t71 * t70;
t75 = pkin(4) * t68 + pkin(6) + qJ(2);
t57 = g(1) * t69 - g(2) * t71;
t60 = pkin(4) * t70 + pkin(3);
t66 = -qJ(5) - pkin(7);
t74 = t60 * t62 - t61 * t66;
t55 = -t62 * t78 + t79;
t53 = t62 * t80 + t77;
t65 = cos(pkin(8));
t72 = pkin(2) * t65 + pkin(1) + t74;
t56 = t62 * t77 + t80;
t54 = -t62 * t79 + t78;
t1 = [(-g(1) * (-t69 * pkin(1) + t71 * qJ(2)) - g(2) * (t71 * pkin(1) + t69 * qJ(2))) * MDP(6) + ((-g(1) * t75 - g(2) * t72) * t71 + (g(1) * t72 - g(2) * t75) * t69) * MDP(24) + (MDP(3) - MDP(5)) * t58 + t89 * (-g(1) * t54 - g(2) * t56) + t88 * (-g(1) * t53 - g(2) * t55) + (MDP(12) * t62 + MDP(4) * t65 - t90 * t61 + MDP(2)) * t57; (-MDP(24) - MDP(6)) * t57; (-g(3) * t74 + t58 * (t60 * t61 + t62 * t66)) * MDP(24) + t90 * (t58 * t62 + t82) + (-t68 * t88 + t70 * t89 + MDP(12)) * t50; t88 * (g(1) * t56 - g(2) * t54 + t70 * t82) + (MDP(24) * pkin(4) + t89) * (-g(1) * t55 + g(2) * t53 + t68 * t82); -t50 * MDP(24);];
taug = t1;
