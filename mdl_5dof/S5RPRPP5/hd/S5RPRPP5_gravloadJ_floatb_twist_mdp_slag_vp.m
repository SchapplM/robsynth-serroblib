% Calculate Gravitation load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:43
% DurationCPUTime: 0.23s
% Computational Cost: add. (94->52), mult. (188->64), div. (0->0), fcn. (147->4), ass. (0->24)
t71 = MDP(12) + MDP(14) + MDP(18);
t70 = MDP(13) - MDP(16) - MDP(19);
t69 = -pkin(1) - pkin(6);
t68 = -pkin(3) - pkin(4);
t57 = cos(qJ(1));
t67 = g(2) * t57;
t54 = sin(qJ(3));
t55 = sin(qJ(1));
t66 = t54 * t55;
t65 = t54 * t57;
t56 = cos(qJ(3));
t64 = t55 * t56;
t61 = qJ(4) * t54;
t63 = pkin(3) * t64 + t55 * t61;
t62 = t57 * pkin(1) + t55 * qJ(2);
t48 = t56 * qJ(4);
t60 = MDP(17) + MDP(21);
t59 = pkin(3) * t66 + t57 * pkin(6) + t62;
t40 = g(1) * t57 + g(2) * t55;
t49 = t57 * qJ(2);
t58 = pkin(3) * t65 - t48 * t57 + t49;
t39 = g(1) * t55 - t67;
t36 = g(1) * t64 - g(3) * t54 - t56 * t67;
t1 = [(-g(1) * (-t55 * pkin(1) + t49) - g(2) * t62) * MDP(6) + (-g(1) * (t55 * t69 + t58) - g(2) * (-t48 * t55 + t59)) * MDP(17) + (-g(1) * (pkin(4) * t65 + t58) - g(2) * (-t57 * qJ(5) + t59) + (-g(1) * (qJ(5) + t69) - g(2) * (t54 * pkin(4) - t48)) * t55) * MDP(21) + (MDP(2) - MDP(4) + MDP(15) - MDP(20)) * t39 + (-t71 * t54 - t70 * t56 + MDP(3) - MDP(5)) * t40; (-MDP(6) - t60) * t39; (-g(1) * t63 - g(3) * (-t54 * pkin(3) + t48) - (-pkin(3) * t56 - t61) * t67) * MDP(17) + (-g(1) * (pkin(4) * t64 + t63) - g(3) * (t54 * t68 + t48) - (t56 * t68 - t61) * t67) * MDP(21) - t71 * t36 + t70 * (g(1) * t66 - g(2) * t65 + g(3) * t56); t60 * t36; t40 * MDP(21);];
taug = t1;
