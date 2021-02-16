% Calculate minimal parameter regressor of potential energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:27
% EndTime: 2021-01-15 14:56:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->25), mult. (84->30), div. (0->0), fcn. (95->6), ass. (0->19)
t42 = cos(qJ(3));
t41 = sin(qJ(3));
t40 = g(3) * qJ(1);
t39 = cos(pkin(7));
t38 = sin(pkin(7));
t37 = t39 * qJ(2);
t24 = -t38 * t41 - t39 * t42;
t25 = -t38 * t42 + t39 * t41;
t36 = g(1) * t25 - g(2) * t24;
t35 = g(1) * t24 + g(2) * t25;
t34 = pkin(1) + pkin(2);
t33 = cos(qJ(4));
t32 = sin(qJ(4));
t31 = qJ(5) + pkin(6);
t30 = t38 * qJ(2);
t29 = t33 * pkin(4) + pkin(3);
t23 = g(3) * t32 + t35 * t33;
t22 = g(3) * t33 - t35 * t32;
t1 = [-t40, -g(1) * (t39 * pkin(1) + t30) - g(2) * (t38 * pkin(1) - t37) - t40, 0, t35, t36, 0, 0, 0, 0, 0, t23, t22, t23, t22, -t36, -g(1) * (-t24 * t29 + t25 * t31 + t34 * t39 + t30) - g(2) * (-t24 * t31 - t25 * t29 + t34 * t38 - t37) - g(3) * (-t32 * pkin(4) - pkin(5) + qJ(1));];
U_reg = t1;
