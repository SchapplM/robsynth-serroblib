% Calculate Gravitation load on the joints for
% S5RPRPR7
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:27
% EndTime: 2021-01-15 12:06:29
% DurationCPUTime: 0.17s
% Computational Cost: add. (136->45), mult. (141->64), div. (0->0), fcn. (119->10), ass. (0->28)
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t60 = MDP(21) * t45 - MDP(22) * t42 + MDP(12);
t39 = qJ(3) + pkin(9);
t35 = sin(t39);
t43 = sin(qJ(3));
t59 = t43 * MDP(11) + MDP(13) * t35;
t58 = g(3) * t35;
t40 = qJ(1) + pkin(8);
t36 = sin(t40);
t57 = t36 * t42;
t56 = t36 * t45;
t38 = cos(t40);
t55 = t38 * t42;
t54 = t38 * t45;
t51 = g(1) * t38 + g(2) * t36;
t50 = g(1) * t36 - g(2) * t38;
t47 = cos(qJ(1));
t46 = cos(qJ(3));
t44 = sin(qJ(1));
t41 = -qJ(4) - pkin(6);
t37 = cos(t39);
t34 = pkin(3) * t46 + pkin(2);
t33 = t37 * t54 + t57;
t32 = -t37 * t55 + t56;
t31 = -t37 * t56 + t55;
t30 = t37 * t57 + t54;
t1 = [(g(1) * t47 + g(2) * t44) * MDP(3) - t51 * MDP(14) + (-g(1) * (-pkin(1) * t44 - t34 * t36 - t38 * t41) - g(2) * (pkin(1) * t47 + t38 * t34 - t36 * t41)) * MDP(15) + (-g(1) * t31 - g(2) * t33) * MDP(21) + (-g(1) * t30 - g(2) * t32) * MDP(22) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t44 - g(2) * t47) + (MDP(10) * t46 + MDP(12) * t37 - t59) * t50; (-MDP(15) - MDP(4)) * g(3); (-t60 * t37 + t59) * g(3) + (t46 * MDP(11) + MDP(13) * t37 + t60 * t35) * t51 + (MDP(15) * pkin(3) + MDP(10)) * (-g(3) * t46 + t51 * t43); -t50 * MDP(15); (-g(1) * t32 + g(2) * t30 + t42 * t58) * MDP(21) + (g(1) * t33 - g(2) * t31 + t45 * t58) * MDP(22);];
taug = t1;
