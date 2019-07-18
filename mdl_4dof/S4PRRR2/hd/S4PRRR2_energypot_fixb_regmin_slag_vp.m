% Calculate minimal parameter regressor of potential energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:24
% EndTime: 2019-07-18 13:27:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->8), mult. (13->13), div. (0->0), fcn. (12->6), ass. (0->9)
t42 = qJ(2) + qJ(3);
t44 = cos(qJ(2));
t43 = sin(qJ(2));
t41 = qJ(4) + t42;
t40 = cos(t42);
t39 = sin(t42);
t38 = cos(t41);
t37 = sin(t41);
t1 = [g(2) * qJ(1), 0, g(1) * t43 - g(3) * t44, g(1) * t44 + g(3) * t43, 0, g(1) * t39 - g(3) * t40, g(1) * t40 + g(3) * t39, 0, g(1) * t37 - g(3) * t38, g(1) * t38 + g(3) * t37;];
U_reg  = t1;
