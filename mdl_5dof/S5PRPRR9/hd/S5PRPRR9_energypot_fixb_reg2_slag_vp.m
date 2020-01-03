% Calculate inertial parameters regressor of potential energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:44
% EndTime: 2019-12-31 17:39:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->36), mult. (97->42), div. (0->0), fcn. (98->8), ass. (0->24)
t42 = pkin(5) + qJ(1);
t57 = g(3) * (-pkin(6) + t42);
t56 = g(3) * t42;
t55 = cos(qJ(4));
t54 = sin(qJ(4));
t53 = pkin(8) + qJ(2);
t36 = cos(t53);
t41 = cos(pkin(8));
t50 = sin(t53);
t52 = t41 * pkin(1) + t36 * pkin(2) + t50 * qJ(3);
t51 = t36 * pkin(3) + t52;
t40 = sin(pkin(8));
t49 = t40 * pkin(1) + t50 * pkin(2) - t36 * qJ(3);
t24 = -t36 * t55 - t50 * t54;
t25 = t36 * t54 - t50 * t55;
t48 = g(1) * t25 - g(2) * t24;
t47 = g(1) * t24 + g(2) * t25;
t46 = -g(1) * t41 - g(2) * t40;
t45 = t50 * pkin(3) + t49;
t44 = cos(qJ(5));
t43 = sin(qJ(5));
t27 = -g(1) * t36 - g(2) * t50;
t26 = g(1) * t50 - g(2) * t36;
t1 = [0, 0, 0, 0, 0, 0, t46, g(1) * t40 - g(2) * t41, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, t27, t26, -g(3), t46 * pkin(1) - t56, 0, 0, 0, 0, 0, 0, t27, -g(3), -t26, -g(1) * t52 - g(2) * t49 - t56, 0, 0, 0, 0, 0, 0, t47, t48, g(3), -g(1) * t51 - g(2) * t45 - t57, 0, 0, 0, 0, 0, 0, g(3) * t43 + t47 * t44, g(3) * t44 - t47 * t43, -t48, -g(1) * (-t24 * pkin(4) + t25 * pkin(7) + t51) - g(2) * (-t25 * pkin(4) - t24 * pkin(7) + t45) - t57;];
U_reg = t1;
