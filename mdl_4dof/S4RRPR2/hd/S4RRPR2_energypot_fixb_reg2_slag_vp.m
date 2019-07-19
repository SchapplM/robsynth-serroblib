% Calculate inertial parameters regressor of potential energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (61->26), mult. (50->31), div. (0->0), fcn. (44->6), ass. (0->16)
t39 = g(3) * (pkin(5) + pkin(4));
t30 = qJ(1) + qJ(2);
t26 = sin(t30);
t27 = cos(t30);
t34 = cos(qJ(1));
t38 = t34 * pkin(1) + t27 * pkin(2) + t26 * qJ(3);
t32 = sin(qJ(1));
t37 = t32 * pkin(1) + t26 * pkin(2) - t27 * qJ(3);
t36 = -g(1) * t34 - g(2) * t32;
t33 = cos(qJ(4));
t31 = sin(qJ(4));
t22 = -g(1) * t27 - g(2) * t26;
t21 = g(1) * t26 - g(2) * t27;
t20 = t26 * t33 - t27 * t31;
t19 = -t26 * t31 - t27 * t33;
t1 = [0, 0, 0, 0, 0, 0, t36, g(1) * t32 - g(2) * t34, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, t22, t21, -g(3), t36 * pkin(1) - t39, 0, 0, 0, 0, 0, 0, t22, -g(3), -t21, -g(1) * t38 - g(2) * t37 - t39, 0, 0, 0, 0, 0, 0, g(1) * t19 - g(2) * t20, -g(1) * t20 - g(2) * t19, g(3), -g(1) * (t27 * pkin(3) + t38) - g(2) * (t26 * pkin(3) + t37) - t39;];
U_reg  = t1;
