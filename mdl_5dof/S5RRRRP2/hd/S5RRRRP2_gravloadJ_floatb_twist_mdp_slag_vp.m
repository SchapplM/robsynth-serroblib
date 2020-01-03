% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:46
% EndTime: 2020-01-03 12:11:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (166->36), mult. (148->50), div. (0->0), fcn. (111->8), ass. (0->21)
t55 = qJ(1) + qJ(2);
t49 = sin(t55);
t51 = cos(t55);
t42 = g(2) * t49 - g(3) * t51;
t54 = qJ(3) + qJ(4);
t48 = sin(t54);
t50 = cos(t54);
t60 = -g(1) * t50 + t42 * t48;
t65 = t60 * MDP(19) + (g(1) * t48 + t42 * t50) * MDP(20);
t58 = cos(qJ(3));
t63 = t58 * pkin(3) + pkin(4) * t50;
t44 = pkin(2) + t63;
t53 = -qJ(5) - pkin(8) - pkin(7);
t64 = t49 * t44 + t51 * t53;
t62 = t51 * t44 - t49 * t53;
t43 = g(2) * t51 + g(3) * t49;
t56 = sin(qJ(3));
t61 = (-MDP(21) + MDP(6)) * t42 + (-t58 * MDP(12) + t56 * MDP(13) - MDP(19) * t50 + MDP(20) * t48 - MDP(5)) * t43;
t59 = cos(qJ(1));
t57 = sin(qJ(1));
t1 = [(-g(2) * t59 - g(3) * t57) * MDP(2) + (g(2) * t57 - g(3) * t59) * MDP(3) + (-g(2) * (t59 * pkin(1) + t62) - g(3) * (t57 * pkin(1) + t64)) * MDP(22) + t61; (-g(2) * t62 - g(3) * t64) * MDP(22) + t61; (-g(1) * t58 + t42 * t56) * MDP(12) + (g(1) * t56 + t42 * t58) * MDP(13) + (-g(1) * t63 - t42 * (-t56 * pkin(3) - pkin(4) * t48)) * MDP(22) + t65; t60 * MDP(22) * pkin(4) + t65; t43 * MDP(22);];
taug = t1;
