% Calculate Gravitation load on the joints for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:48
% EndTime: 2019-03-09 02:33:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (166->47), mult. (182->62), div. (0->0), fcn. (156->10), ass. (0->25)
t54 = pkin(10) + qJ(4);
t49 = qJ(5) + t54;
t46 = cos(t49);
t68 = g(3) * t46;
t57 = sin(qJ(6));
t60 = cos(qJ(1));
t67 = t57 * t60;
t58 = sin(qJ(1));
t66 = t58 * t57;
t59 = cos(qJ(6));
t65 = t58 * t59;
t64 = t59 * t60;
t63 = t60 * pkin(1) + t58 * qJ(2);
t43 = g(1) * t58 - g(2) * t60;
t45 = sin(t49);
t62 = (t43 * t45 + t68) * MDP(24) + (-t59 * MDP(30) + t57 * MDP(31) - MDP(23)) * (-g(3) * t45 + t43 * t46);
t44 = g(1) * t60 + g(2) * t58;
t51 = t60 * qJ(2);
t48 = cos(t54);
t47 = sin(t54);
t42 = t45 * t64 - t66;
t41 = t45 * t67 + t65;
t40 = t45 * t65 + t67;
t39 = -t45 * t66 + t64;
t1 = [(-g(1) * (-t58 * pkin(1) + t51) - g(2) * t63) * MDP(6) + (-g(1) * (t51 + (-pkin(1) - qJ(3)) * t58) - g(2) * (qJ(3) * t60 + t63)) * MDP(10) + (-g(1) * t42 - g(2) * t40) * MDP(30) + (g(1) * t41 - g(2) * t39) * MDP(31) + (MDP(2) - MDP(4) + MDP(9)) * t43 + (-t47 * MDP(16) - t48 * MDP(17) - MDP(23) * t45 - MDP(24) * t46 - MDP(7) * sin(pkin(10)) - MDP(8) * cos(pkin(10)) + MDP(3) - MDP(5)) * t44; (-MDP(10) - MDP(6)) * t43; -t44 * MDP(10); (g(3) * t47 - t43 * t48) * MDP(16) + (g(3) * t48 + t43 * t47) * MDP(17) + t62; t62; (-g(1) * t39 - g(2) * t41 + t57 * t68) * MDP(30) + (g(1) * t40 - g(2) * t42 + t59 * t68) * MDP(31);];
taug  = t1;
