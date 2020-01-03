% Calculate Gravitation load on the joints for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (144->42), mult. (128->52), div. (0->0), fcn. (95->8), ass. (0->21)
t42 = qJ(1) + pkin(7);
t37 = sin(t42);
t39 = cos(t42);
t53 = g(1) * t39 + g(2) * t37;
t59 = MDP(14) + MDP(16);
t58 = MDP(15) - MDP(18);
t46 = sin(qJ(1));
t55 = t46 * pkin(1);
t54 = MDP(19) + MDP(8);
t52 = g(1) * t37 - g(2) * t39;
t41 = pkin(8) + qJ(4);
t36 = sin(t41);
t38 = cos(t41);
t50 = t38 * pkin(4) + t36 * qJ(5);
t44 = cos(pkin(8));
t48 = t44 * pkin(3) + pkin(2) + t50;
t47 = cos(qJ(1));
t45 = -pkin(6) - qJ(3);
t40 = t47 * pkin(1);
t29 = -g(3) * t38 + t36 * t53;
t1 = [(g(1) * t47 + g(2) * t46) * MDP(3) + (-g(1) * (-t37 * pkin(2) + t39 * qJ(3) - t55) - g(2) * (t39 * pkin(2) + t37 * qJ(3) + t40)) * MDP(8) + (g(1) * t55 - g(2) * t40 + (g(1) * t45 - g(2) * t48) * t39 + (g(1) * t48 + g(2) * t45) * t37) * MDP(19) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t46 - g(2) * t47) - (MDP(7) + MDP(17)) * t53 + (MDP(5) * t44 - MDP(6) * sin(pkin(8)) - t58 * t36 + t59 * t38) * t52; (-MDP(4) - t54) * g(3); -t54 * t52; (-g(3) * t50 + t53 * (pkin(4) * t36 - qJ(5) * t38)) * MDP(19) + t58 * (g(3) * t36 + t38 * t53) + t59 * t29; -t29 * MDP(19);];
taug = t1;
