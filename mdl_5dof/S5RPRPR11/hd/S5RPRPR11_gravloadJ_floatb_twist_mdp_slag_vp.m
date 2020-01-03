% Calculate Gravitation load on the joints for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:56
% DurationCPUTime: 0.21s
% Computational Cost: add. (153->44), mult. (194->62), div. (0->0), fcn. (177->8), ass. (0->24)
t54 = pkin(8) + qJ(3);
t51 = sin(t54);
t52 = cos(t54);
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t75 = -t51 * t60 + t52 * t58;
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t49 = g(1) * t61 + g(2) * t59;
t74 = MDP(13) + MDP(15);
t73 = MDP(14) - MDP(17);
t39 = t75 * t59;
t65 = t51 * t58 + t52 * t60;
t40 = t65 * t59;
t41 = t75 * t61;
t42 = t65 * t61;
t72 = (g(1) * t41 + g(2) * t39 + g(3) * t65) * MDP(24) + (g(1) * t42 + g(2) * t40 - g(3) * t75) * MDP(25);
t48 = g(1) * t59 - g(2) * t61;
t67 = t52 * pkin(3) + t51 * qJ(4);
t56 = cos(pkin(8));
t64 = t56 * pkin(2) + pkin(1) + t67;
t57 = -pkin(6) - qJ(2);
t37 = -g(3) * t52 + t49 * t51;
t1 = [(-g(1) * (-t59 * pkin(1) + t61 * qJ(2)) - g(2) * (t61 * pkin(1) + t59 * qJ(2))) * MDP(7) + ((g(1) * t57 - g(2) * t64) * t61 + (g(1) * t64 + g(2) * t57) * t59) * MDP(18) + (g(1) * t40 - g(2) * t42) * MDP(24) + (-g(1) * t39 + g(2) * t41) * MDP(25) + (MDP(3) - MDP(6) - MDP(16)) * t49 + (MDP(4) * t56 - MDP(5) * sin(pkin(8)) - t73 * t51 + t74 * t52 + MDP(2)) * t48; (-MDP(18) - MDP(7)) * t48; (-g(3) * t67 + t49 * (pkin(3) * t51 - qJ(4) * t52)) * MDP(18) + t73 * (g(3) * t51 + t49 * t52) + t74 * t37 - t72; -t37 * MDP(18); t72;];
taug = t1;
