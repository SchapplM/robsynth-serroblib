% Calculate Gravitation load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:27
% EndTime: 2019-07-18 17:22:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (104->37), mult. (162->56), div. (0->0), fcn. (141->8), ass. (0->23)
t45 = cos(qJ(2));
t57 = pkin(1) * t45;
t40 = qJ(2) + qJ(4);
t38 = sin(t40);
t56 = g(3) * t38;
t41 = sin(qJ(5));
t43 = sin(qJ(1));
t54 = t43 * t41;
t44 = cos(qJ(5));
t53 = t43 * t44;
t46 = cos(qJ(1));
t52 = t46 * t41;
t51 = t46 * t44;
t39 = cos(t40);
t49 = g(1) * t46 + g(2) * t43;
t50 = (t49 * t39 + t56) * MDP(19) + (MDP(25) * t44 - MDP(26) * t41 + MDP(18)) * (-g(3) * t39 + t49 * t38);
t36 = g(1) * t43 - g(2) * t46;
t42 = sin(qJ(2));
t35 = t39 * t51 + t54;
t34 = -t39 * t52 + t53;
t33 = -t39 * t53 + t52;
t32 = t39 * t54 + t51;
t1 = [(-g(1) * (t46 * qJ(3) - t43 * t57) - g(2) * (t43 * qJ(3) + t46 * t57)) * MDP(12) + (-g(1) * t33 - g(2) * t35) * MDP(25) + (-g(1) * t32 - g(2) * t34) * MDP(26) + (MDP(3) - MDP(11)) * t49 + (-t42 * MDP(10) + MDP(18) * t39 - MDP(19) * t38 + t45 * MDP(9) + MDP(2)) * t36; (g(3) * t42 + t49 * t45) * MDP(10) + t50 + (MDP(12) * pkin(1) + MDP(9)) * (-g(3) * t45 + t49 * t42); -t36 * MDP(12); t50; (-g(1) * t34 + g(2) * t32 + t41 * t56) * MDP(25) + (g(1) * t35 - g(2) * t33 + t44 * t56) * MDP(26);];
taug  = t1;
