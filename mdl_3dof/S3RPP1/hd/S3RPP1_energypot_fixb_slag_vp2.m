% Calculate potential energy for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3RPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energypot_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energypot_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_energypot_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPP1_energypot_fixb_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:32
% EndTime: 2018-11-14 10:13:32
% DurationCPUTime: 0.05s
% Computational Cost: add. (31->21), mult. (41->15), div. (0->0), fcn. (22->2), ass. (0->6)
t21 = m(3) + m(4);
t19 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2);
t18 = -m(4) * qJ(3) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t1 = (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - t21) * pkin(3)) * g(3) + (-mrSges(1,2) + (t21 * qJ(2) - t19) * t17 + (-t21 * pkin(1) + t18) * t16) * g(2) + (t19 * t16 + t18 * t17 - mrSges(1,1) - t21 * (t17 * pkin(1) + t16 * qJ(2))) * g(1);
U  = t1;
