% Calculate inertial parameters regressor of potential energy for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:40
% EndTime: 2018-11-14 13:51:40
% DurationCPUTime: 0.04s
% Computational Cost: add. (64->27), mult. (40->28), div. (0->0), fcn. (30->6), ass. (0->16)
t39 = pkin(5) + pkin(4);
t38 = g(3) * (qJ(3) + t39);
t32 = qJ(1) + qJ(2);
t27 = sin(t32);
t33 = sin(qJ(1));
t37 = t33 * pkin(1) + pkin(2) * t27;
t28 = cos(t32);
t34 = cos(qJ(1));
t36 = t34 * pkin(1) + pkin(2) * t28;
t35 = -g(1) * t34 - g(2) * t33;
t26 = pkin(6) + t32;
t23 = cos(t26);
t22 = sin(t26);
t21 = -g(1) * t23 - g(2) * t22;
t20 = g(1) * t22 - g(2) * t23;
t1 = [0, 0, 0, 0, 0, 0, t35, g(1) * t33 - g(2) * t34, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t27, g(1) * t27 - g(2) * t28, -g(3), t35 * pkin(1) - g(3) * t39, 0, 0, 0, 0, 0, 0, t21, t20, -g(3), -g(1) * t36 - g(2) * t37 - t38, 0, 0, 0, 0, 0, 0, t21, -g(3), -t20, -g(1) * (t23 * pkin(3) + t22 * qJ(4) + t36) - g(2) * (t22 * pkin(3) - t23 * qJ(4) + t37) - t38;];
U_reg  = t1;
