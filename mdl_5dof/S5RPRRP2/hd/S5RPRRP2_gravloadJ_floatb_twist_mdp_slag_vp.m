% Calculate Gravitation load on the joints for
% S5RPRRP2
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
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:30
% DurationCPUTime: 0.08s
% Computational Cost: add. (128->30), mult. (97->39), div. (0->0), fcn. (67->8), ass. (0->16)
t37 = qJ(1) + pkin(8);
t36 = qJ(3) + t37;
t33 = sin(t36);
t34 = cos(t36);
t41 = cos(qJ(4));
t35 = t41 * pkin(4) + pkin(3);
t38 = -qJ(5) - pkin(7);
t47 = t33 * t35 + t34 * t38;
t46 = -t33 * t38 + t34 * t35;
t28 = g(2) * t33 - g(3) * t34;
t29 = g(2) * t34 + g(3) * t33;
t39 = sin(qJ(4));
t45 = (-MDP(15) + MDP(7)) * t28 + (-t41 * MDP(13) + t39 * MDP(14) - MDP(6)) * t29;
t42 = cos(qJ(1));
t40 = sin(qJ(1));
t1 = [(g(2) * t40 - g(3) * t42) * MDP(3) + (-g(2) * (pkin(2) * cos(t37) + t42 * pkin(1) + t46) - g(3) * (pkin(2) * sin(t37) + t40 * pkin(1) + t47)) * MDP(16) + t45 + (pkin(1) * MDP(4) + MDP(2)) * (-g(2) * t42 - g(3) * t40); (-MDP(16) - MDP(4)) * g(1); (-g(2) * t46 - g(3) * t47) * MDP(16) + t45; (g(1) * t39 + t28 * t41) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (-g(1) * t41 + t28 * t39); t29 * MDP(16);];
taug = t1;
