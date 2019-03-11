% Calculate potential energy for
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
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S3RRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energypot_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (43->27), mult. (37->20), div. (0->0), fcn. (18->6), ass. (0->12)
t24 = pkin(4) + pkin(3);
t19 = qJ(1) + qJ(2);
t23 = -m(4) * pkin(2) - mrSges(3,1);
t22 = -mrSges(2,1) + (-m(3) - m(4)) * pkin(1);
t21 = cos(qJ(1));
t20 = sin(qJ(1));
t18 = qJ(3) + t19;
t17 = cos(t19);
t16 = sin(t19);
t15 = cos(t18);
t14 = sin(t18);
t1 = (-mrSges(1,3) - m(2) * pkin(3) - mrSges(2,3) - m(3) * t24 - mrSges(3,3) - m(4) * (pkin(5) + t24) - mrSges(4,3)) * g(3) + (-t14 * mrSges(4,1) - t21 * mrSges(2,2) - t17 * mrSges(3,2) - t15 * mrSges(4,2) + t23 * t16 + t22 * t20 - mrSges(1,2)) * g(2) + (-t15 * mrSges(4,1) + t20 * mrSges(2,2) + t16 * mrSges(3,2) + t14 * mrSges(4,2) + t23 * t17 + t22 * t21 - mrSges(1,1)) * g(1);
U  = t1;
