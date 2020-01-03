% Calculate Gravitation load on the joints for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:21
% EndTime: 2019-12-31 18:57:21
% DurationCPUTime: 0.22s
% Computational Cost: add. (86->47), mult. (178->65), div. (0->0), fcn. (153->6), ass. (0->25)
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t77 = -g(1) * t56 + g(2) * t59;
t76 = MDP(13) - MDP(21);
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t40 = -g(3) * t55 - t58 * t77;
t69 = g(3) * t58;
t54 = sin(qJ(4));
t68 = t56 * t54;
t57 = cos(qJ(4));
t67 = t56 * t57;
t66 = t59 * t54;
t65 = t59 * t57;
t63 = pkin(4) * t54 + pkin(6);
t62 = g(2) * (t59 * pkin(1) + t56 * qJ(2));
t48 = t57 * pkin(4) + pkin(3);
t53 = -qJ(5) - pkin(7);
t60 = t55 * t48 + t58 * t53;
t43 = t55 * t66 + t67;
t41 = -t55 * t68 + t65;
t50 = t59 * qJ(2);
t44 = t55 * t65 - t68;
t42 = t55 * t67 + t66;
t1 = [(-g(1) * (-t56 * pkin(1) + t50) - t62) * MDP(6) + (-g(1) * t44 - g(2) * t42) * MDP(19) + (g(1) * t43 - g(2) * t41) * MDP(20) + (-g(1) * t50 - t62 + (-g(1) * t60 - g(2) * t63) * t59 + (-g(1) * (-pkin(1) - t63) - g(2) * t60) * t56) * MDP(22) - (MDP(2) - MDP(4)) * t77 + (-t55 * MDP(12) - t76 * t58 + MDP(3) - MDP(5)) * (g(1) * t59 + g(2) * t56); -(-MDP(22) - MDP(6)) * t77; (g(3) * t60 + t77 * (t48 * t58 - t53 * t55)) * MDP(22) + t76 * (-t55 * t77 + t69) + (-MDP(19) * t57 + MDP(20) * t54 - MDP(12)) * t40; (g(1) * t42 - g(2) * t44 + t57 * t69) * MDP(20) + (pkin(4) * MDP(22) + MDP(19)) * (-g(1) * t41 - g(2) * t43 + t54 * t69); t40 * MDP(22);];
taug = t1;
