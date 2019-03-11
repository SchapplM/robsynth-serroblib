% Calculate inertial parameters regressor of potential energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% U_reg [1x(3*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3RPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energypot_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->19), mult. (37->24), div. (0->0), fcn. (34->4), ass. (0->12)
t24 = g(3) * pkin(3);
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t23 = t21 * pkin(1) + t19 * qJ(2);
t22 = t19 * pkin(1) - t21 * qJ(2);
t20 = cos(qJ(3));
t18 = sin(qJ(3));
t14 = -g(1) * t21 - g(2) * t19;
t13 = g(1) * t19 - g(2) * t21;
t12 = -t21 * t18 + t19 * t20;
t11 = -t19 * t18 - t21 * t20;
t1 = [0, 0, 0, 0, 0, 0, t14, t13, -g(3), -t24, 0, 0, 0, 0, 0, 0, t14, -g(3), -t13, -g(1) * t23 - g(2) * t22 - t24, 0, 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t12, -g(1) * t12 - g(2) * t11, g(3), -g(1) * (t21 * pkin(2) + t23) - g(2) * (t19 * pkin(2) + t22) - g(3) * (-pkin(4) + pkin(3));];
U_reg  = t1;
