% Calculate Gravitation load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:02
% DurationCPUTime: 0.16s
% Computational Cost: add. (105->40), mult. (128->58), div. (0->0), fcn. (111->8), ass. (0->23)
t44 = qJ(1) + pkin(7);
t41 = sin(t44);
t59 = g(1) * t41;
t58 = MDP(12) + MDP(8);
t42 = cos(t44);
t50 = cos(qJ(1));
t57 = t50 * pkin(1) + t42 * pkin(2) + t41 * qJ(3);
t48 = sin(qJ(1));
t56 = -t48 * pkin(1) + t42 * qJ(3);
t37 = -g(1) * t42 - g(2) * t41;
t55 = -g(2) * t42 + t59;
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t53 = pkin(3) * t46 + qJ(4) * t45;
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t52 = t45 * t49 - t46 * t47;
t51 = t45 * t47 + t46 * t49;
t33 = t51 * t42;
t32 = t52 * t42;
t31 = t51 * t41;
t30 = t52 * t41;
t1 = [(g(1) * t50 + g(2) * t48) * MDP(3) + (-g(1) * (-t41 * pkin(2) + t56) - g(2) * t57) * MDP(8) + (-g(1) * t56 - g(2) * (t42 * t53 + t57) - (-pkin(2) - t53) * t59) * MDP(12) + (g(1) * t31 - g(2) * t33) * MDP(18) + (g(1) * t30 - g(2) * t32) * MDP(19) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t48 - g(2) * t50) + (MDP(7) + MDP(10)) * t37 + ((MDP(5) + MDP(9)) * t46 + (-MDP(6) + MDP(11)) * t45) * t55; (-MDP(4) - t58) * g(3); -t58 * t55; (g(3) * t46 + t37 * t45) * MDP(12); (-g(1) * t32 - g(2) * t30 + g(3) * t51) * MDP(18) + (g(1) * t33 + g(2) * t31 + g(3) * t52) * MDP(19);];
taug = t1;
