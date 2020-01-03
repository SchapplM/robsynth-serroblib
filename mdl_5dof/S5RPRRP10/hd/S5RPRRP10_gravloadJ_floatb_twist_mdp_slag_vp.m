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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:06
% DurationCPUTime: 0.26s
% Computational Cost: add. (137->49), mult. (185->67), div. (0->0), fcn. (159->8), ass. (0->27)
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t48 = g(1) * t62 + g(2) * t60;
t79 = MDP(14) - MDP(22);
t54 = pkin(8) + qJ(3);
t51 = sin(t54);
t52 = cos(t54);
t40 = -g(3) * t52 + t48 * t51;
t73 = g(3) * t51;
t59 = sin(qJ(4));
t71 = t59 * t62;
t70 = t60 * t59;
t61 = cos(qJ(4));
t69 = t60 * t61;
t68 = t61 * t62;
t66 = pkin(4) * t59 + pkin(6) + qJ(2);
t47 = g(1) * t60 - g(2) * t62;
t50 = pkin(4) * t61 + pkin(3);
t57 = -qJ(5) - pkin(7);
t65 = t50 * t52 - t51 * t57;
t45 = -t52 * t71 + t69;
t43 = t52 * t70 + t68;
t56 = cos(pkin(8));
t63 = pkin(2) * t56 + pkin(1) + t65;
t46 = t52 * t68 + t70;
t44 = -t52 * t69 + t71;
t1 = [(-g(1) * (-t60 * pkin(1) + qJ(2) * t62) - g(2) * (pkin(1) * t62 + t60 * qJ(2))) * MDP(7) + (-g(1) * t44 - g(2) * t46) * MDP(20) + (-g(1) * t43 - g(2) * t45) * MDP(21) + ((-g(1) * t66 - g(2) * t63) * t62 + (g(1) * t63 - g(2) * t66) * t60) * MDP(23) + (MDP(3) - MDP(6)) * t48 + (t52 * MDP(13) + MDP(4) * t56 - MDP(5) * sin(pkin(8)) - t79 * t51 + MDP(2)) * t47; (-MDP(23) - MDP(7)) * t47; (-g(3) * t65 + t48 * (t50 * t51 + t52 * t57)) * MDP(23) + t79 * (t48 * t52 + t73) + (MDP(20) * t61 - MDP(21) * t59 + MDP(13)) * t40; (g(1) * t46 - g(2) * t44 + t61 * t73) * MDP(21) + (pkin(4) * MDP(23) + MDP(20)) * (-g(1) * t45 + g(2) * t43 + t59 * t73); -t40 * MDP(23);];
taug = t1;
