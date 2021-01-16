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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:53
% EndTime: 2021-01-15 20:08:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (170->38), mult. (138->49), div. (0->0), fcn. (99->8), ass. (0->25)
t69 = -MDP(13) - MDP(15);
t68 = MDP(14) + MDP(16);
t55 = qJ(1) + qJ(2);
t50 = pkin(8) + t55;
t45 = sin(t50);
t46 = cos(t50);
t51 = sin(t55);
t47 = pkin(2) * t51;
t59 = cos(qJ(4));
t49 = t59 * pkin(4) + pkin(3);
t56 = -qJ(5) - pkin(7);
t67 = t45 * t49 + t46 * t56 + t47;
t52 = cos(t55);
t48 = pkin(2) * t52;
t66 = -t45 * t56 + t46 * t49 + t48;
t65 = g(2) * t46 + g(3) * t45;
t64 = g(2) * t45 - g(3) * t46;
t63 = -g(2) * t52 - g(3) * t51;
t57 = sin(qJ(4));
t62 = -t64 * MDP(17) + (g(2) * t51 - g(3) * t52) * MDP(6) + t63 * MDP(5) + (t68 * t57 + t69 * t59) * t65;
t60 = cos(qJ(1));
t58 = sin(qJ(1));
t54 = t60 * pkin(1);
t53 = t58 * pkin(1);
t1 = [(-g(2) * t60 - g(3) * t58) * MDP(2) + (g(2) * t58 - g(3) * t60) * MDP(3) + (-g(2) * (t48 + t54) - g(3) * (t47 + t53)) * MDP(7) + (-g(2) * (t54 + t66) - g(3) * (t53 + t67)) * MDP(18) + t62; (-g(2) * t66 - g(3) * t67) * MDP(18) + t63 * MDP(7) * pkin(2) + t62; (-MDP(18) - MDP(7)) * g(1); t68 * (g(1) * t57 + t64 * t59) + (MDP(18) * pkin(4) - t69) * (-g(1) * t59 + t64 * t57); t65 * MDP(18);];
taug = t1;
