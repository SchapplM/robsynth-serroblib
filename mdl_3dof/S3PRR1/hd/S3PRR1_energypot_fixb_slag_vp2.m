% Calculate potential energy for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S3PRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_energypot_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_energypot_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:49
% EndTime: 2019-03-08 18:03:49
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->26), mult. (29->18), div. (0->0), fcn. (10->4), ass. (0->6)
t14 = cos(qJ(2));
t13 = sin(qJ(2));
t12 = qJ(2) + qJ(3);
t11 = cos(t12);
t10 = sin(t12);
t1 = (-mrSges(1,3) + mrSges(2,2) - m(3) * pkin(3) - mrSges(3,3) - m(4) * (pkin(4) + pkin(3)) - mrSges(4,3)) * g(3) + (-t10 * mrSges(4,1) - t14 * mrSges(3,2) - t11 * mrSges(4,2) - mrSges(1,2) - mrSges(2,3) + (-m(4) * pkin(2) - mrSges(3,1)) * t13 + (-m(2) - m(3) - m(4)) * qJ(1)) * g(2) + (-mrSges(1,1) - mrSges(2,1) - m(3) * pkin(1) - t14 * mrSges(3,1) + t13 * mrSges(3,2) - m(4) * (t14 * pkin(2) + pkin(1)) - t11 * mrSges(4,1) + t10 * mrSges(4,2)) * g(1);
U  = t1;
