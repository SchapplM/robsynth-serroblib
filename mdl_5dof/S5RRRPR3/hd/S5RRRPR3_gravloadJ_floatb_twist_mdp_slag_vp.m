% Calculate Gravitation load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:29
% EndTime: 2020-01-03 12:09:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (160->32), mult. (132->43), div. (0->0), fcn. (99->8), ass. (0->19)
t48 = qJ(1) + qJ(2);
t46 = sin(t48);
t47 = cos(t48);
t37 = g(2) * t46 - g(3) * t47;
t45 = qJ(3) + pkin(9) + qJ(5);
t42 = sin(t45);
t43 = cos(t45);
t58 = (-g(1) * t43 + t37 * t42) * MDP(21) + (g(1) * t42 + t37 * t43) * MDP(22);
t52 = cos(qJ(3));
t44 = t52 * pkin(3) + pkin(2);
t49 = -qJ(4) - pkin(7);
t57 = t46 * t44 + t47 * t49;
t56 = t47 * t44 - t46 * t49;
t38 = g(2) * t47 + g(3) * t46;
t50 = sin(qJ(3));
t55 = (-MDP(14) + MDP(6)) * t37 + (-t52 * MDP(12) + t50 * MDP(13) - t43 * MDP(21) + t42 * MDP(22) - MDP(5)) * t38;
t53 = cos(qJ(1));
t51 = sin(qJ(1));
t1 = [(-g(2) * t53 - g(3) * t51) * MDP(2) + (g(2) * t51 - g(3) * t53) * MDP(3) + (-g(2) * (t53 * pkin(1) + t56) - g(3) * (t51 * pkin(1) + t57)) * MDP(15) + t55; (-g(2) * t56 - g(3) * t57) * MDP(15) + t55; (g(1) * t50 + t37 * t52) * MDP(13) + t58 + (MDP(15) * pkin(3) + MDP(12)) * (-g(1) * t52 + t37 * t50); t38 * MDP(15); t58;];
taug = t1;
