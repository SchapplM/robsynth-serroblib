% Calculate Gravitation load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (71->33), mult. (166->51), div. (0->0), fcn. (157->6), ass. (0->20)
t42 = sin(qJ(4));
t43 = sin(qJ(2));
t45 = cos(qJ(4));
t46 = cos(qJ(2));
t51 = t46 * t42 - t43 * t45;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t40 = g(1) * t47 + g(2) * t44;
t61 = MDP(9) + MDP(11);
t60 = MDP(10) - MDP(13);
t32 = t51 * t44;
t38 = t43 * t42 + t46 * t45;
t33 = t38 * t44;
t34 = t51 * t47;
t35 = t38 * t47;
t59 = (g(1) * t34 + g(2) * t32 + g(3) * t38) * MDP(20) + (g(1) * t35 + g(2) * t33 - g(3) * t51) * MDP(21);
t53 = t46 * pkin(2) + t43 * qJ(3);
t50 = pkin(1) + t53;
t30 = -g(3) * t46 + t40 * t43;
t1 = [((-g(1) * pkin(5) - g(2) * t50) * t47 + (-g(2) * pkin(5) + g(1) * t50) * t44) * MDP(14) + (g(1) * t33 - g(2) * t35) * MDP(20) + (-g(1) * t32 + g(2) * t34) * MDP(21) + (MDP(3) - MDP(12)) * t40 + (-t60 * t43 + t61 * t46 + MDP(2)) * (g(1) * t44 - g(2) * t47); (-g(3) * t53 + t40 * (pkin(2) * t43 - qJ(3) * t46)) * MDP(14) + t60 * (g(3) * t43 + t40 * t46) + t61 * t30 - t59; -t30 * MDP(14); t59;];
taug = t1;
