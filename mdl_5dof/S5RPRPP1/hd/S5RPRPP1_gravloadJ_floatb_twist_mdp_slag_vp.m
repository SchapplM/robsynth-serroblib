% Calculate Gravitation load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:14
% EndTime: 2019-12-31 18:09:15
% DurationCPUTime: 0.20s
% Computational Cost: add. (140->46), mult. (134->58), div. (0->0), fcn. (99->8), ass. (0->21)
t42 = qJ(3) + pkin(8);
t36 = sin(t42);
t38 = cos(t42);
t62 = t38 * pkin(4) + t36 * qJ(5);
t43 = qJ(1) + pkin(7);
t37 = sin(t43);
t39 = cos(t43);
t55 = g(1) * t39 + g(2) * t37;
t47 = cos(qJ(3));
t40 = t47 * pkin(3);
t35 = t40 + pkin(2);
t48 = cos(qJ(1));
t58 = t48 * pkin(1) + t39 * t35;
t56 = MDP(13) + MDP(17);
t54 = g(1) * t37 - g(2) * t39;
t44 = -qJ(4) - pkin(6);
t46 = sin(qJ(1));
t52 = -t46 * pkin(1) - t39 * t44;
t45 = sin(qJ(3));
t31 = -g(3) * t38 + t55 * t36;
t1 = [(g(1) * t48 + g(2) * t46) * MDP(3) + (-g(1) * (-t37 * t35 + t52) - g(2) * (-t37 * t44 + t58)) * MDP(13) + (-g(1) * t52 - g(2) * (t62 * t39 + t58) + (-g(1) * (-t35 - t62) + g(2) * t44) * t37) * MDP(17) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t46 - g(2) * t48) - (MDP(12) + MDP(15)) * t55 + (t47 * MDP(10) - t45 * MDP(11) + t38 * MDP(14) + t36 * MDP(16)) * t54; (-MDP(4) - t56) * g(3); (g(3) * t45 + t55 * t47) * MDP(11) + t31 * MDP(14) + (-g(3) * t36 - t55 * t38) * MDP(16) + (-g(3) * (t40 + t62) + t55 * (pkin(3) * t45 + pkin(4) * t36 - qJ(5) * t38)) * MDP(17) + (pkin(3) * MDP(13) + MDP(10)) * (-g(3) * t47 + t55 * t45); -t56 * t54; -t31 * MDP(17);];
taug = t1;
