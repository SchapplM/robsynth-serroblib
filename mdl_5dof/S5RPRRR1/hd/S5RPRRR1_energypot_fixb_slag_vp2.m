% Calculate potential energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:33
% EndTime: 2019-12-05 18:08:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (62->35), mult. (125->32), div. (0->0), fcn. (115->8), ass. (0->16)
t44 = -mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + (m(3) + m(4) + m(5) + m(6)) * qJ(2);
t31 = sin(qJ(1));
t34 = cos(qJ(3));
t43 = t31 * t34;
t35 = cos(qJ(1));
t42 = t34 * t35;
t41 = mrSges(5,2) - mrSges(6,3);
t28 = sin(qJ(5));
t32 = cos(qJ(5));
t38 = -mrSges(6,1) * t32 + mrSges(6,2) * t28 - mrSges(5,1);
t37 = t28 * mrSges(6,1) + t32 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t30 = sin(qJ(3));
t36 = -t34 * mrSges(4,1) - t37 * t30 - mrSges(2,1) - mrSges(3,1);
t33 = cos(qJ(4));
t29 = sin(qJ(4));
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t37 * t34 + (t41 * t29 + t38 * t33 - mrSges(4,1)) * t30) * g(3) + (-mrSges(1,2) + t41 * (t29 * t43 + t33 * t35) + t38 * (-t29 * t35 + t33 * t43) + t36 * t31 + t44 * t35) * g(2) + (-mrSges(1,1) + t41 * (t29 * t42 - t31 * t33) + t38 * (t31 * t29 + t33 * t42) + t36 * t35 - t44 * t31) * g(1);
U = t1;
