% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x8]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:49
% EndTime: 2019-03-08 18:21:49
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->11), mult. (14->14), div. (0->0), fcn. (10->4), ass. (0->6)
t56 = cos(qJ(2));
t55 = sin(qJ(2));
t54 = qJ(2) + pkin(6) + qJ(4);
t53 = cos(t54);
t52 = sin(t54);
t1 = [-g(2) * qJ(1), 0, -g(1) * t56 - g(2) * t55, g(1) * t55 - g(2) * t56, -g(1) * (t56 * pkin(2) + pkin(1)) - g(2) * (t55 * pkin(2) + qJ(1)) - g(3) * (qJ(3) + pkin(4)) 0, -g(1) * t53 - g(2) * t52, g(1) * t52 - g(2) * t53;];
U_reg  = t1;
