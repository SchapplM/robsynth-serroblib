% Calculate Gravitation load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:12:55
% DurationCPUTime: 0.30s
% Computational Cost: add. (165->49), mult. (198->58), div. (0->0), fcn. (156->6), ass. (0->27)
t58 = cos(pkin(7));
t56 = pkin(7) + qJ(3);
t53 = sin(t56);
t54 = cos(t56);
t67 = t54 * pkin(3) + t53 * qJ(4);
t77 = pkin(2) * t58 + pkin(1) + t67;
t61 = cos(qJ(1));
t76 = g(2) * t61;
t75 = MDP(13) - MDP(16) + MDP(21);
t74 = MDP(14) - MDP(17) - MDP(20);
t60 = sin(qJ(1));
t71 = g(1) * t61;
t47 = g(2) * t60 + t71;
t73 = t47 * t53;
t72 = pkin(3) * t53;
t59 = -pkin(6) - qJ(2);
t69 = pkin(4) - t59;
t66 = qJ(4) * t54;
t65 = qJ(5) * t54;
t64 = -MDP(18) - MDP(22);
t63 = t77 * t76;
t46 = g(1) * t60 - t76;
t44 = t61 * t66;
t42 = t60 * t66;
t39 = g(3) * t53 + t47 * t54;
t38 = -g(3) * t54 + t73;
t1 = [(-g(1) * (-t60 * pkin(1) + qJ(2) * t61) - g(2) * (pkin(1) * t61 + t60 * qJ(2))) * MDP(7) + (t59 * t71 - t63 + (g(1) * t77 + g(2) * t59) * t60) * MDP(18) + (-t63 + (-g(1) * t69 - g(2) * t65) * t61 + (-g(1) * (-t77 - t65) - g(2) * t69) * t60) * MDP(22) + (MDP(3) - MDP(6) - MDP(15) - MDP(19)) * t47 + (t58 * MDP(4) - MDP(5) * sin(pkin(7)) + MDP(2) + t75 * t54 - t74 * t53) * t46; (-MDP(7) + t64) * t46; (-g(1) * (-t61 * t72 + t44) - g(2) * (-t60 * t72 + t42) - g(3) * t67) * MDP(18) + (-g(1) * t44 - g(2) * t42 - g(3) * (t65 + t67) + (pkin(3) + qJ(5)) * t73) * MDP(22) + t74 * t39 + t75 * t38; t64 * t38; -t39 * MDP(22);];
taug = t1;
