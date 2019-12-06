% Calculate Gravitation load on the joints for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:07
% EndTime: 2019-12-05 15:36:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (124->40), mult. (186->59), div. (0->0), fcn. (171->8), ass. (0->25)
t67 = MDP(11) + MDP(13);
t66 = MDP(12) - MDP(15);
t46 = sin(pkin(7));
t47 = cos(pkin(7));
t57 = g(1) * t47 + g(2) * t46;
t45 = qJ(2) + pkin(8);
t43 = sin(t45);
t63 = g(3) * t43;
t48 = sin(qJ(4));
t62 = t46 * t48;
t50 = cos(qJ(4));
t61 = t46 * t50;
t60 = t47 * t48;
t59 = t47 * t50;
t58 = MDP(16) + MDP(5);
t56 = pkin(4) * t50 + qJ(5) * t48 + pkin(3);
t44 = cos(t45);
t36 = t44 * t62 + t59;
t38 = t44 * t60 - t61;
t33 = g(1) * t38 + g(2) * t36 + t48 * t63;
t51 = cos(qJ(2));
t49 = sin(qJ(2));
t39 = t44 * t59 + t62;
t37 = t44 * t61 - t60;
t1 = [(-MDP(1) - t58) * g(3); (g(3) * t49 + t51 * t57) * MDP(4) + (-t44 * t57 - t63) * MDP(14) + (-g(3) * (t51 * pkin(2) + t43 * pkin(6) + t44 * t56) + t57 * (pkin(2) * t49 - pkin(6) * t44 + t43 * t56)) * MDP(16) + (pkin(2) * MDP(5) + MDP(3)) * (-g(3) * t51 + t49 * t57) + (-t66 * t48 + t67 * t50) * (-g(3) * t44 + t43 * t57); t58 * (-g(1) * t46 + g(2) * t47); (-g(1) * (-t38 * pkin(4) + t39 * qJ(5)) - g(2) * (-t36 * pkin(4) + t37 * qJ(5)) - (-pkin(4) * t48 + qJ(5) * t50) * t63) * MDP(16) + t66 * (g(1) * t39 + g(2) * t37 + t50 * t63) + t67 * t33; -t33 * MDP(16);];
taug = t1;
