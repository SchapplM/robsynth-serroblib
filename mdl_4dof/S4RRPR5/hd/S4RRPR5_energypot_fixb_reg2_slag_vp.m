% Calculate inertial parameters regressor of potential energy for
% S4RRPR5
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
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:32
% EndTime: 2019-12-31 17:03:32
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->27), mult. (50->28), div. (0->0), fcn. (40->6), ass. (0->15)
t41 = pkin(5) + pkin(4);
t45 = g(3) * t41;
t36 = qJ(1) + qJ(2);
t32 = sin(t36);
t33 = cos(t36);
t40 = cos(qJ(1));
t44 = t40 * pkin(1) + t33 * pkin(2) + t32 * qJ(3);
t38 = sin(qJ(1));
t43 = t38 * pkin(1) + t32 * pkin(2) - t33 * qJ(3);
t27 = g(1) * t32 - g(2) * t33;
t42 = -g(1) * t40 - g(2) * t38;
t39 = cos(qJ(4));
t37 = sin(qJ(4));
t28 = g(1) * t33 + g(2) * t32;
t1 = [0, 0, 0, 0, 0, 0, t42, g(1) * t38 - g(2) * t40, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t28, t27, -g(3), t42 * pkin(1) - t45, 0, 0, 0, 0, 0, 0, -g(3), t28, -t27, -g(1) * t44 - g(2) * t43 - t45, 0, 0, 0, 0, 0, 0, -g(3) * t39 - t27 * t37, g(3) * t37 - t27 * t39, -t28, -g(1) * (t33 * pkin(6) + t44) - g(2) * (t32 * pkin(6) + t43) - g(3) * (pkin(3) + t41);];
U_reg = t1;
