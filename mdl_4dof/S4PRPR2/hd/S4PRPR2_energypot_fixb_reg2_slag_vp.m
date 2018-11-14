% Calculate inertial parameters regressor of potential energy for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:21
% EndTime: 2018-11-14 14:02:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (41->26), mult. (28->25), div. (0->0), fcn. (18->6), ass. (0->13)
t29 = cos(qJ(2));
t33 = t29 * pkin(2) + pkin(1);
t32 = g(2) * qJ(1);
t31 = qJ(3) + pkin(4);
t28 = sin(qJ(2));
t30 = t28 * pkin(2) + qJ(1);
t27 = qJ(2) + pkin(6);
t24 = qJ(4) + t27;
t23 = cos(t27);
t22 = sin(t27);
t21 = cos(t24);
t20 = sin(t24);
t1 = [0, 0, 0, 0, 0, 0, -g(1), g(3), -g(2), -t32, 0, 0, 0, 0, 0, 0, -g(1) * t29 - g(2) * t28, g(1) * t28 - g(2) * t29, -g(3), -g(1) * pkin(1) - g(3) * pkin(4) - t32, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t22, g(1) * t22 - g(2) * t23, -g(3), -g(1) * t33 - g(2) * t30 - g(3) * t31, 0, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t20, g(1) * t20 - g(2) * t21, -g(3), -g(1) * (pkin(3) * t23 + t33) - g(2) * (pkin(3) * t22 + t30) - g(3) * (pkin(5) + t31);];
U_reg  = t1;
