% Calculate Gravitation load on the joints for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:38
% EndTime: 2019-12-31 18:24:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (120->42), mult. (155->64), div. (0->0), fcn. (133->8), ass. (0->25)
t44 = qJ(1) + pkin(8);
t41 = sin(t44);
t42 = cos(t44);
t58 = g(1) * t42 + g(2) * t41;
t67 = MDP(10) - MDP(13);
t66 = MDP(11) - MDP(14);
t49 = cos(qJ(3));
t61 = g(3) * t49;
t45 = sin(qJ(5));
t46 = sin(qJ(3));
t60 = t45 * t46;
t48 = cos(qJ(5));
t59 = t46 * t48;
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t56 = g(1) * t47 - g(2) * t50;
t55 = t49 * pkin(3) + t46 * qJ(4);
t53 = pkin(2) + t55;
t52 = t56 * pkin(1);
t38 = -t41 * t60 + t42 * t48;
t37 = t41 * t59 + t42 * t45;
t36 = t41 * t48 + t42 * t60;
t35 = -t41 * t45 + t42 * t59;
t33 = t58 * t46 - t61;
t1 = [t56 * MDP(2) + (g(1) * t50 + g(2) * t47) * MDP(3) + MDP(4) * t52 - t58 * MDP(12) + (t52 + (-g(1) * pkin(6) - g(2) * t53) * t42 + (-g(2) * pkin(6) + g(1) * t53) * t41) * MDP(15) + (-g(1) * t38 - g(2) * t36) * MDP(21) + (g(1) * t37 - g(2) * t35) * MDP(22) + (-t66 * t46 + t67 * t49) * (g(1) * t41 - g(2) * t42); (-MDP(15) - MDP(4)) * g(3); (-g(3) * t55 + t58 * (pkin(3) * t46 - qJ(4) * t49)) * MDP(15) + t67 * t33 + (-MDP(21) * t45 - MDP(22) * t48 + t66) * (g(3) * t46 + t58 * t49); -t33 * MDP(15); (-g(1) * t35 - g(2) * t37 + t48 * t61) * MDP(21) + (g(1) * t36 - g(2) * t38 - t45 * t61) * MDP(22);];
taug = t1;
