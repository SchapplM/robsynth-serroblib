% Calculate minimal parameter regressor of potential energy for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:13:59
% EndTime: 2021-01-15 15:13:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (76->33), mult. (93->47), div. (0->0), fcn. (94->8), ass. (0->26)
t53 = cos(qJ(4));
t41 = t53 * pkin(4) + pkin(3);
t46 = qJ(2) + pkin(8);
t43 = sin(t46);
t44 = cos(t46);
t49 = -qJ(5) - pkin(6);
t65 = t41 * t44 - t43 * t49;
t64 = g(3) * t43;
t47 = sin(pkin(7));
t51 = sin(qJ(4));
t61 = t47 * t51;
t60 = t47 * t53;
t48 = cos(pkin(7));
t59 = t48 * t51;
t58 = t48 * t53;
t54 = cos(qJ(2));
t42 = t54 * pkin(2) + pkin(1);
t50 = -qJ(3) - pkin(5);
t57 = t47 * t42 + t48 * t50;
t52 = sin(qJ(2));
t56 = t52 * pkin(2) + qJ(1);
t55 = g(1) * t48 + g(2) * t47;
t39 = t48 * t42;
t37 = -g(1) * (t44 * t58 + t61) - g(2) * (t44 * t60 - t59) - t53 * t64;
t36 = -g(1) * (-t44 * t59 + t60) - g(2) * (-t44 * t61 - t58) + t51 * t64;
t1 = [-g(3) * qJ(1), 0, -g(3) * t52 - t55 * t54, -g(3) * t54 + t55 * t52, -g(1) * (-t47 * t50 + t39) - g(2) * t57 - g(3) * t56, 0, 0, 0, 0, 0, t37, t36, t37, t36, g(3) * t44 - t55 * t43, -g(1) * (t65 * t48 + t39) - g(2) * (-pkin(4) * t59 + t57) - g(3) * (t43 * t41 + t44 * t49 + t56) + (-g(1) * (pkin(4) * t51 - t50) - g(2) * t65) * t47;];
U_reg = t1;
