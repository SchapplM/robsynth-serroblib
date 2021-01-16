% Calculate minimal parameter regressor of potential energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:39
% EndTime: 2021-01-15 16:14:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->24), mult. (45->26), div. (0->0), fcn. (41->8), ass. (0->15)
t35 = pkin(8) + qJ(2);
t34 = qJ(3) + t35;
t29 = sin(t34);
t30 = cos(t34);
t39 = g(1) * t30 + g(2) * t29;
t38 = cos(qJ(4));
t37 = sin(qJ(4));
t36 = qJ(5) + pkin(7);
t33 = cos(t35);
t32 = sin(t35);
t31 = t38 * pkin(4) + pkin(3);
t28 = g(1) * t29 - g(2) * t30;
t27 = -g(3) * t37 - t39 * t38;
t26 = -g(3) * t38 + t39 * t37;
t1 = [-g(3) * qJ(1), 0, -g(1) * t33 - g(2) * t32, g(1) * t32 - g(2) * t33, 0, -t39, t28, 0, 0, 0, 0, 0, t27, t26, t27, t26, -t28, -g(1) * (t30 * t31 + t36 * t29 + cos(pkin(8)) * pkin(1) + pkin(2) * t33) - g(2) * (t29 * t31 - t30 * t36 + pkin(2) * t32 + sin(pkin(8)) * pkin(1)) - g(3) * (t37 * pkin(4) + pkin(5) + pkin(6) + qJ(1));];
U_reg = t1;
