% Calculate Gravitation load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:05
% EndTime: 2019-12-31 16:54:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (72->32), mult. (109->48), div. (0->0), fcn. (98->8), ass. (0->21)
t31 = pkin(7) + qJ(3);
t29 = sin(t31);
t47 = t29 * MDP(14);
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t46 = t36 * MDP(20) - t34 * MDP(21) + MDP(13);
t45 = g(3) * t29;
t35 = sin(qJ(1));
t44 = t35 * t34;
t43 = t35 * t36;
t37 = cos(qJ(1));
t42 = t37 * t34;
t41 = t37 * t36;
t38 = g(1) * t37 + g(2) * t35;
t27 = g(1) * t35 - g(2) * t37;
t30 = cos(t31);
t26 = t30 * t41 + t44;
t25 = -t30 * t42 + t43;
t24 = -t30 * t43 + t42;
t23 = t30 * t44 + t41;
t1 = [(-g(1) * (-t35 * pkin(1) + t37 * qJ(2)) - g(2) * (t37 * pkin(1) + t35 * qJ(2))) * MDP(7) + (-g(1) * t24 - g(2) * t26) * MDP(20) + (-g(1) * t23 - g(2) * t25) * MDP(21) + (MDP(3) - MDP(6)) * t38 + (t30 * MDP(13) - t47 + MDP(4) * cos(pkin(7)) - MDP(5) * sin(pkin(7)) + MDP(2)) * t27; -t27 * MDP(7); (-t46 * t30 + t47) * g(3) + (MDP(14) * t30 + t46 * t29) * t38; (-g(1) * t25 + g(2) * t23 + t34 * t45) * MDP(20) + (g(1) * t26 - g(2) * t24 + t36 * t45) * MDP(21);];
taug = t1;
