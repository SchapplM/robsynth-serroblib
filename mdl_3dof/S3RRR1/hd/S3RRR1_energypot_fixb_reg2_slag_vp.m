% Calculate inertial parameters regressor of potential energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% U_reg [1x(3*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S3RRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energypot_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:59
% EndTime: 2018-11-14 10:15:59
% DurationCPUTime: 0.03s
% Computational Cost: add. (33->19), mult. (25->22), div. (0->0), fcn. (18->6), ass. (0->11)
t24 = pkin(4) + pkin(3);
t20 = qJ(1) + qJ(2);
t21 = sin(qJ(1));
t22 = cos(qJ(1));
t23 = -g(1) * t22 - g(2) * t21;
t19 = qJ(3) + t20;
t18 = cos(t20);
t17 = sin(t20);
t16 = cos(t19);
t15 = sin(t19);
t1 = [0, 0, 0, 0, 0, 0, t23, g(1) * t21 - g(2) * t22, -g(3), -g(3) * pkin(3), 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t17, g(1) * t17 - g(2) * t18, -g(3), t23 * pkin(1) - g(3) * t24, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t15, g(1) * t15 - g(2) * t16, -g(3), -g(1) * (t22 * pkin(1) + pkin(2) * t18) - g(2) * (t21 * pkin(1) + pkin(2) * t17) - g(3) * (pkin(5) + t24);];
U_reg  = t1;
