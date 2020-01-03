% Calculate Gravitation load on the joints for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:17
% DurationCPUTime: 0.15s
% Computational Cost: add. (78->43), mult. (133->60), div. (0->0), fcn. (113->8), ass. (0->25)
t46 = sin(qJ(3));
t59 = pkin(3) * t46;
t43 = qJ(3) + pkin(8);
t38 = cos(t43);
t57 = g(3) * t38;
t45 = sin(qJ(5));
t47 = sin(qJ(1));
t56 = t47 * t45;
t48 = cos(qJ(5));
t55 = t47 * t48;
t50 = cos(qJ(1));
t54 = t50 * t45;
t53 = t50 * t48;
t52 = t50 * pkin(1) + t47 * qJ(2);
t36 = g(1) * t50 + g(2) * t47;
t35 = g(1) * t47 - g(2) * t50;
t49 = cos(qJ(3));
t44 = -qJ(4) - pkin(6);
t40 = t50 * qJ(2);
t37 = sin(t43);
t34 = t37 * t53 - t56;
t33 = t37 * t54 + t55;
t32 = t37 * t55 + t54;
t31 = -t37 * t56 + t53;
t1 = [(-g(1) * (-t47 * pkin(1) + t40) - g(2) * t52) * MDP(6) + (-g(1) * (t50 * t59 + t40 + (-pkin(1) + t44) * t47) - g(2) * (-t50 * t44 + t47 * t59 + t52)) * MDP(15) + (-g(1) * t34 - g(2) * t32) * MDP(21) + (g(1) * t33 - g(2) * t31) * MDP(22) + (MDP(2) - MDP(4) + MDP(14)) * t35 + (-t46 * MDP(12) - t49 * MDP(13) + MDP(3) - MDP(5)) * t36; (-MDP(15) - MDP(6)) * t35; (g(3) * t49 + t35 * t46) * MDP(13) + (MDP(15) * pkin(3) + MDP(12)) * (g(3) * t46 - t35 * t49) + (-MDP(21) * t48 + MDP(22) * t45) * (-g(3) * t37 + t35 * t38); -t36 * MDP(15); (-g(1) * t31 - g(2) * t33 + t45 * t57) * MDP(21) + (g(1) * t32 - g(2) * t34 + t48 * t57) * MDP(22);];
taug = t1;
