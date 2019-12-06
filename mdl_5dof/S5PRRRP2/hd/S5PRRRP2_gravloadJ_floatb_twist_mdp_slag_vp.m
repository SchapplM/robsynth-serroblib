% Calculate Gravitation load on the joints for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:58
% EndTime: 2019-12-05 16:41:59
% DurationCPUTime: 0.09s
% Computational Cost: add. (190->33), mult. (142->41), div. (0->0), fcn. (109->6), ass. (0->20)
t55 = sin(qJ(4));
t56 = cos(qJ(4));
t60 = t56 * pkin(4) + t55 * qJ(5);
t70 = -pkin(3) - t60;
t69 = MDP(13) + MDP(15);
t68 = MDP(14) - MDP(17);
t54 = pkin(8) + qJ(2);
t53 = qJ(3) + t54;
t49 = sin(t53);
t50 = cos(t53);
t43 = g(1) * t50 + g(2) * t49;
t67 = g(1) * t49;
t62 = t49 * pkin(7) - t70 * t50;
t58 = (-MDP(16) + MDP(7)) * t43 + (-t68 * t55 + t69 * t56 + MDP(6)) * (-g(2) * t50 + t67);
t57 = t70 * t67;
t52 = cos(t54);
t51 = sin(t54);
t47 = t50 * pkin(7);
t32 = -g(3) * t56 + t43 * t55;
t1 = [(-MDP(1) - MDP(18)) * g(3); (g(1) * t51 - g(2) * t52) * MDP(3) + (g(1) * t52 + g(2) * t51) * MDP(4) + (-g(1) * (-pkin(2) * t51 + t47) - g(2) * (pkin(2) * t52 + t62) - t57) * MDP(18) + t58; (-g(1) * t47 - g(2) * t62 - t57) * MDP(18) + t58; (-g(3) * t60 + t43 * (pkin(4) * t55 - qJ(5) * t56)) * MDP(18) + t68 * (g(3) * t55 + t43 * t56) + t69 * t32; -t32 * MDP(18);];
taug = t1;
