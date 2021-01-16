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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:01
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:01:02
% EndTime: 2021-01-16 00:01:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (222->38), mult. (192->50), div. (0->0), fcn. (147->8), ass. (0->23)
t72 = MDP(19) + MDP(21);
t71 = MDP(20) + MDP(22);
t60 = qJ(3) + qJ(4);
t56 = cos(t60);
t64 = cos(qJ(3));
t70 = t64 * pkin(3) + pkin(4) * t56;
t50 = -pkin(2) - t70;
t61 = qJ(1) + qJ(2);
t55 = sin(t61);
t57 = cos(t61);
t59 = -qJ(5) - pkin(8) - pkin(7);
t69 = -t57 * t50 - t55 * t59;
t68 = -t50 * t55 + t57 * t59;
t48 = g(2) * t55 - g(3) * t57;
t54 = sin(t60);
t35 = -g(1) * t56 + t48 * t54;
t67 = t72 * t35 + t71 * (g(1) * t54 + t48 * t56);
t49 = g(2) * t57 + g(3) * t55;
t62 = sin(qJ(3));
t66 = (-MDP(23) + MDP(6)) * t48 + (-t64 * MDP(12) + t62 * MDP(13) + t71 * t54 - t72 * t56 - MDP(5)) * t49;
t65 = cos(qJ(1));
t63 = sin(qJ(1));
t1 = [(-g(2) * t65 - g(3) * t63) * MDP(2) + (g(2) * t63 - g(3) * t65) * MDP(3) + (-g(2) * (pkin(1) * t65 + t69) - g(3) * (pkin(1) * t63 + t68)) * MDP(24) + t66; (-g(2) * t69 - g(3) * t68) * MDP(24) + t66; (-g(1) * t64 + t48 * t62) * MDP(12) + (g(1) * t62 + t48 * t64) * MDP(13) + (-g(1) * t70 + t48 * (pkin(3) * t62 + pkin(4) * t54)) * MDP(24) + t67; t35 * MDP(24) * pkin(4) + t67; t49 * MDP(24);];
taug = t1;
