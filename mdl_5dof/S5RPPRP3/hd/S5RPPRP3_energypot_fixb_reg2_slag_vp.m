% Calculate inertial parameters regressor of potential energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:01
% EndTime: 2019-12-31 17:51:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (93->37), mult. (78->37), div. (0->0), fcn. (65->6), ass. (0->21)
t49 = sin(qJ(4));
t59 = pkin(4) * t49;
t48 = qJ(2) + pkin(5);
t58 = g(3) * t48;
t46 = qJ(1) + pkin(7);
t42 = sin(t46);
t50 = sin(qJ(1));
t57 = t50 * pkin(1) + t42 * pkin(2);
t56 = pkin(3) + t48;
t43 = cos(t46);
t52 = cos(qJ(1));
t55 = t52 * pkin(1) + t43 * pkin(2) + t42 * qJ(3);
t54 = -t43 * qJ(3) + t57;
t37 = g(1) * t42 - g(2) * t43;
t53 = -g(1) * t52 - g(2) * t50;
t51 = cos(qJ(4));
t47 = -qJ(5) - pkin(6);
t38 = g(1) * t43 + g(2) * t42;
t36 = g(3) * t49 - t37 * t51;
t35 = -g(3) * t51 - t37 * t49;
t1 = [0, 0, 0, 0, 0, 0, t53, g(1) * t50 - g(2) * t52, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t38, t37, -g(3), t53 * pkin(1) - t58, 0, 0, 0, 0, 0, 0, -g(3), t38, -t37, -g(1) * t55 - g(2) * t54 - t58, 0, 0, 0, 0, 0, 0, t35, t36, -t38, -g(1) * (t43 * pkin(6) + t55) - g(2) * (t42 * pkin(6) + t54) - g(3) * t56, 0, 0, 0, 0, 0, 0, t35, t36, -t38, -g(1) * (t42 * t59 - t43 * t47 + t55) - g(2) * (-t42 * t47 + (-qJ(3) - t59) * t43 + t57) - g(3) * (t51 * pkin(4) + t56);];
U_reg = t1;
