% Calculate potential energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:48
% EndTime: 2019-12-05 17:02:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (69->36), mult. (92->27), div. (0->0), fcn. (67->8), ass. (0->17)
t42 = -m(3) - m(4);
t41 = m(5) + m(6);
t34 = cos(qJ(3));
t40 = pkin(2) * t34;
t39 = mrSges(5,2) - mrSges(6,3);
t30 = sin(qJ(5));
t33 = cos(qJ(5));
t38 = mrSges(6,1) * t33 - mrSges(6,2) * t30 + mrSges(5,1);
t37 = t30 * mrSges(6,1) + t33 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t29 = qJ(3) + qJ(4);
t27 = sin(t29);
t28 = cos(t29);
t31 = sin(qJ(3));
t36 = -mrSges(4,1) * t34 + mrSges(4,2) * t31 + t39 * t27 - t38 * t28 - mrSges(3,1);
t35 = cos(qJ(2));
t32 = sin(qJ(2));
t1 = (-mrSges(1,3) - mrSges(2,3) - t41 * (t32 * t40 + qJ(1)) + (-m(2) + t42) * qJ(1) + t37 * t35 + t36 * t32) * g(3) + (t34 * mrSges(4,2) - mrSges(1,2) - mrSges(2,2) + mrSges(3,3) + t39 * t28 + (t41 * pkin(2) + mrSges(4,1)) * t31 + t38 * t27) * g(2) + (-mrSges(1,1) - mrSges(2,1) - t41 * (t35 * t40 + pkin(1)) + t42 * pkin(1) - t37 * t32 + t36 * t35) * g(1);
U = t1;
