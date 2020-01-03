% Calculate Gravitation load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:02
% EndTime: 2020-01-03 11:34:02
% DurationCPUTime: 0.07s
% Computational Cost: add. (151->30), mult. (104->39), div. (0->0), fcn. (74->10), ass. (0->15)
t43 = qJ(1) + pkin(8);
t41 = qJ(3) + t43;
t37 = sin(t41);
t38 = cos(t41);
t51 = t38 * pkin(3) + t37 * qJ(4);
t50 = t37 * pkin(3) - t38 * qJ(4);
t33 = g(2) * t38 + g(3) * t37;
t32 = g(2) * t37 - g(3) * t38;
t42 = pkin(9) + qJ(5);
t39 = sin(t42);
t40 = cos(t42);
t48 = (-MDP(10) + MDP(7)) * t32 + (-t40 * MDP(17) + t39 * MDP(18) - MDP(8) * cos(pkin(9)) + MDP(9) * sin(pkin(9)) - MDP(6)) * t33;
t47 = cos(qJ(1));
t46 = sin(qJ(1));
t1 = [(g(2) * t46 - g(3) * t47) * MDP(3) + (-g(2) * (pkin(2) * cos(t43) + t47 * pkin(1) + t51) - g(3) * (pkin(2) * sin(t43) + t46 * pkin(1) + t50)) * MDP(11) + t48 + (pkin(1) * MDP(4) + MDP(2)) * (-g(2) * t47 - g(3) * t46); (-MDP(11) - MDP(4)) * g(1); (-g(2) * t51 - g(3) * t50) * MDP(11) + t48; t33 * MDP(11); (-g(1) * t40 + t32 * t39) * MDP(17) + (g(1) * t39 + t32 * t40) * MDP(18);];
taug = t1;
