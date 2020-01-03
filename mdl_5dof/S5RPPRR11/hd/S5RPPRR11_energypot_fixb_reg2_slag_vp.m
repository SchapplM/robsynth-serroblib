% Calculate inertial parameters regressor of potential energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (62->42), mult. (98->48), div. (0->0), fcn. (89->6), ass. (0->25)
t63 = g(3) * pkin(5);
t62 = pkin(2) + pkin(5);
t47 = cos(qJ(4));
t61 = g(3) * t47;
t43 = sin(qJ(5));
t45 = sin(qJ(1));
t60 = t45 * t43;
t46 = cos(qJ(5));
t59 = t45 * t46;
t48 = cos(qJ(1));
t58 = t48 * t43;
t57 = t48 * t46;
t56 = t48 * pkin(1) + t45 * qJ(2);
t55 = pkin(3) + t62;
t54 = t48 * qJ(3) + t56;
t53 = t45 * pkin(1) - t48 * qJ(2);
t52 = t45 * qJ(3) + t53;
t44 = sin(qJ(4));
t51 = pkin(4) * t44 - pkin(7) * t47;
t35 = g(1) * t48 + g(2) * t45;
t50 = -t45 * pkin(6) + t54;
t49 = t48 * pkin(6) + t52;
t34 = g(1) * t45 - g(2) * t48;
t33 = -g(3) * t44 + t35 * t47;
t1 = [0, 0, 0, 0, 0, 0, -t35, t34, -g(3), -t63, 0, 0, 0, 0, 0, 0, -g(3), t35, -t34, -g(1) * t56 - g(2) * t53 - t63, 0, 0, 0, 0, 0, 0, -g(3), -t34, -t35, -g(1) * t54 - g(2) * t52 - g(3) * t62, 0, 0, 0, 0, 0, 0, -t35 * t44 - t61, -t33, t34, -g(1) * t50 - g(2) * t49 - g(3) * t55, 0, 0, 0, 0, 0, 0, -g(1) * (t44 * t57 - t60) - g(2) * (t44 * t59 + t58) - t46 * t61, -g(1) * (-t44 * t58 - t59) - g(2) * (-t44 * t60 + t57) + t43 * t61, t33, -g(1) * (t51 * t48 + t50) - g(2) * (t51 * t45 + t49) - g(3) * (t47 * pkin(4) + t44 * pkin(7) + t55);];
U_reg = t1;
