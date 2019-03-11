% Calculate minimal parameter regressor of potential energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3RRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:52
% EndTime: 2019-03-08 18:06:52
% DurationCPUTime: 0.02s
% Computational Cost: add. (25->13), mult. (21->17), div. (0->0), fcn. (18->4), ass. (0->8)
t36 = cos(qJ(1));
t35 = sin(qJ(1));
t34 = qJ(1) + qJ(2);
t33 = cos(t34);
t32 = sin(t34);
t31 = -g(1) * t33 - g(2) * t32;
t30 = g(1) * t32 - g(2) * t33;
t1 = [0, -g(1) * t36 - g(2) * t35, g(1) * t35 - g(2) * t36, 0, t31, t30, t31, -t30, -g(1) * (t36 * pkin(1) + t33 * pkin(2) + t32 * qJ(3)) - g(2) * (t35 * pkin(1) + t32 * pkin(2) - t33 * qJ(3)) - g(3) * (pkin(4) + pkin(3));];
U_reg  = t1;
