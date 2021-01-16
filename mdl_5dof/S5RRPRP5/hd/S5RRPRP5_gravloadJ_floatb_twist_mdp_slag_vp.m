% Calculate Gravitation load on the joints for
% S5RRPRP5
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:19
% EndTime: 2021-01-15 20:18:21
% DurationCPUTime: 0.28s
% Computational Cost: add. (224->56), mult. (209->69), div. (0->0), fcn. (164->8), ass. (0->29)
t63 = qJ(2) + pkin(8);
t58 = cos(t63);
t67 = cos(qJ(2));
t60 = t67 * pkin(2);
t59 = qJ(4) + t63;
t54 = sin(t59);
t55 = cos(t59);
t75 = t55 * pkin(4) + t54 * qJ(5);
t79 = pkin(3) * t58 + t60 + t75;
t78 = MDP(20) + MDP(22);
t77 = MDP(21) - MDP(24);
t76 = pkin(4) * t54;
t64 = -qJ(3) - pkin(6);
t73 = qJ(5) * t55;
t57 = sin(t63);
t65 = sin(qJ(2));
t72 = -pkin(2) * t65 - pkin(3) * t57 - t76;
t66 = sin(qJ(1));
t68 = cos(qJ(1));
t50 = g(1) * t68 + g(2) * t66;
t41 = -g(3) * t55 + t50 * t54;
t71 = t77 * (g(3) * t54 + t50 * t55) + t78 * t41;
t49 = g(1) * t66 - g(2) * t68;
t70 = pkin(1) + t79;
t62 = -pkin(7) + t64;
t56 = t60 + pkin(1);
t48 = t68 * t73;
t47 = t66 * t73;
t1 = [(-g(1) * (-t66 * t56 - t64 * t68) - g(2) * (t56 * t68 - t66 * t64)) * MDP(14) + ((g(1) * t62 - g(2) * t70) * t68 + (g(1) * t70 + g(2) * t62) * t66) * MDP(25) + (MDP(3) - MDP(13) - MDP(23)) * t50 + (-t65 * MDP(10) + t58 * MDP(11) - t57 * MDP(12) + t67 * MDP(9) - t77 * t54 + t78 * t55 + MDP(2)) * t49; (g(3) * t65 + t50 * t67) * MDP(10) + (-g(3) * t58 + t50 * t57) * MDP(11) + (g(3) * t57 + t50 * t58) * MDP(12) + (-g(1) * (t72 * t68 + t48) - g(2) * (t72 * t66 + t47) - g(3) * t79) * MDP(25) + t71 + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t67 + t50 * t65); (-MDP(14) - MDP(25)) * t49; (-g(1) * (-t68 * t76 + t48) - g(2) * (-t66 * t76 + t47) - g(3) * t75) * MDP(25) + t71; -t41 * MDP(25);];
taug = t1;
