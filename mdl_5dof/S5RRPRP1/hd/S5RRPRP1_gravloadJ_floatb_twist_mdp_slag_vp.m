% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:12
% EndTime: 2020-01-03 11:59:12
% DurationCPUTime: 0.07s
% Computational Cost: add. (132->36), mult. (106->49), div. (0->0), fcn. (73->8), ass. (0->23)
t48 = qJ(1) + qJ(2);
t43 = pkin(8) + t48;
t38 = sin(t43);
t39 = cos(t43);
t44 = sin(t48);
t40 = pkin(2) * t44;
t52 = cos(qJ(4));
t42 = t52 * pkin(4) + pkin(3);
t49 = -qJ(5) - pkin(7);
t60 = t38 * t42 + t39 * t49 + t40;
t45 = cos(t48);
t41 = pkin(2) * t45;
t59 = -t38 * t49 + t39 * t42 + t41;
t50 = sin(qJ(4));
t55 = -g(2) * t45 - g(3) * t44;
t56 = g(2) * t38 - g(3) * t39;
t57 = g(2) * t39 + g(3) * t38;
t58 = -t56 * MDP(15) + (g(2) * t44 - g(3) * t45) * MDP(6) + t55 * MDP(5) + (-t52 * MDP(13) + t50 * MDP(14)) * t57;
t53 = cos(qJ(1));
t51 = sin(qJ(1));
t47 = t53 * pkin(1);
t46 = t51 * pkin(1);
t1 = [(-g(2) * t53 - g(3) * t51) * MDP(2) + (g(2) * t51 - g(3) * t53) * MDP(3) + (-g(2) * (t41 + t47) - g(3) * (t40 + t46)) * MDP(7) + (-g(2) * (t47 + t59) - g(3) * (t46 + t60)) * MDP(16) + t58; (-g(2) * t59 - g(3) * t60) * MDP(16) + t55 * MDP(7) * pkin(2) + t58; (-MDP(16) - MDP(7)) * g(1); (g(1) * t50 + t56 * t52) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (-g(1) * t52 + t56 * t50); t57 * MDP(16);];
taug = t1;
