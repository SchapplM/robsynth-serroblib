% Calculate Gravitation load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:17
% DurationCPUTime: 0.15s
% Computational Cost: add. (113->40), mult. (283->63), div. (0->0), fcn. (315->8), ass. (0->30)
t86 = MDP(11) + MDP(13);
t85 = MDP(12) - MDP(15);
t60 = sin(pkin(8));
t82 = g(3) * t60;
t64 = sin(qJ(4));
t81 = t60 * t64;
t66 = cos(qJ(4));
t80 = t60 * t66;
t67 = cos(qJ(3));
t79 = t60 * t67;
t61 = sin(pkin(7));
t65 = sin(qJ(3));
t78 = t61 * t65;
t77 = t61 * t67;
t63 = cos(pkin(7));
t76 = t63 * t65;
t75 = t63 * t67;
t74 = MDP(16) + MDP(2);
t62 = cos(pkin(8));
t50 = t62 * t77 - t76;
t52 = t62 * t75 + t78;
t73 = -g(1) * t52 - g(2) * t50;
t45 = t50 * t64 - t61 * t80;
t47 = t52 * t64 - t63 * t80;
t53 = t62 * t66 + t64 * t79;
t71 = g(1) * t47 + g(2) * t45 + g(3) * t53;
t54 = -t62 * t64 + t66 * t79;
t48 = t52 * t66 + t63 * t81;
t46 = t50 * t66 + t61 * t81;
t1 = [(-MDP(1) - t74) * g(3); t74 * (-g(1) * t61 + g(2) * t63); (-t67 * t82 + t73) * MDP(16) * pkin(6) + (MDP(5) - MDP(14)) * (g(3) * t79 - t73) + (MDP(4) + (pkin(4) * t66 + qJ(5) * t64 + pkin(3)) * MDP(16) + t86 * t66 - t85 * t64) * (t65 * t82 - g(2) * (-t62 * t78 - t75) - g(1) * (-t62 * t76 + t77)); (-g(1) * (-t47 * pkin(4) + t48 * qJ(5)) - g(2) * (-t45 * pkin(4) + t46 * qJ(5)) - g(3) * (-t53 * pkin(4) + t54 * qJ(5))) * MDP(16) + t86 * t71 + t85 * (g(1) * t48 + g(2) * t46 + g(3) * t54); -t71 * MDP(16);];
taug = t1;
