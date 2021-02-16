% Calculate Gravitation load on the joints for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:31
% EndTime: 2021-01-15 19:36:32
% DurationCPUTime: 0.20s
% Computational Cost: add. (169->51), mult. (218->66), div. (0->0), fcn. (197->10), ass. (0->28)
t62 = qJ(2) + pkin(8);
t58 = sin(t62);
t59 = cos(t62);
t73 = sin(qJ(5));
t74 = cos(qJ(5));
t46 = -t58 * t73 - t59 * t74;
t47 = t58 * t74 - t59 * t73;
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t71 = g(1) * t69 + g(2) * t67;
t82 = -(g(3) * t46 + t71 * t47) * MDP(24) - (-g(3) * t47 + t71 * t46) * MDP(25);
t81 = MDP(11) + MDP(15);
t80 = MDP(12) - MDP(17);
t75 = g(3) * t59;
t65 = -qJ(3) - pkin(6);
t72 = t67 * t65;
t51 = g(1) * t67 - g(2) * t69;
t68 = cos(qJ(2));
t66 = sin(qJ(2));
t64 = cos(pkin(8));
t63 = sin(pkin(8));
t60 = t68 * pkin(2);
t57 = t60 + pkin(1);
t56 = t65 * t69;
t50 = -pkin(3) * t63 + qJ(4) * t64;
t49 = pkin(3) * t64 + qJ(4) * t63 + pkin(2);
t40 = t49 * t68 + t50 * t66 + pkin(1);
t1 = [(-g(1) * (-t57 * t67 - t56) - g(2) * (t57 * t69 - t72)) * MDP(14) + (-g(1) * (-t40 * t67 - t56) - g(2) * (t40 * t69 - t72)) * MDP(18) + (MDP(3) - MDP(13) - MDP(16)) * t71 + (-t66 * MDP(10) - t46 * MDP(24) + t47 * MDP(25) + t68 * MDP(9) - t80 * t58 + t81 * t59 + MDP(2)) * t51; (g(3) * t66 + t71 * t68) * MDP(10) + (-g(3) * (pkin(3) * t59 + qJ(4) * t58 + t60) - t71 * (-t49 * t66 + t50 * t68)) * MDP(18) + (pkin(2) * MDP(14) + MDP(9)) * (-g(3) * t68 + t71 * t66) + t80 * (g(3) * t58 + t71 * t59) + t81 * (t71 * t58 - t75) - t82; (-MDP(14) - MDP(18)) * t51; (t75 - t71 * (t63 * t68 + t64 * t66)) * MDP(18); t82;];
taug = t1;
