% Calculate potential energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:34
% EndTime: 2019-12-31 19:52:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (114->44), mult. (101->34), div. (0->0), fcn. (70->6), ass. (0->18)
t60 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t59 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t58 = -m(5) - m(6);
t57 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t56 = t60 * t43 + t59 * t45 - mrSges(3,2) + mrSges(4,3);
t47 = pkin(6) + pkin(5);
t44 = sin(qJ(1));
t40 = t44 * pkin(1);
t46 = cos(qJ(1));
t41 = t46 * pkin(1);
t42 = qJ(1) + qJ(2);
t38 = sin(t42);
t55 = t38 * pkin(2) + t40;
t39 = cos(t42);
t52 = t39 * pkin(2) + t38 * qJ(3) + t41;
t1 = (-m(2) * pkin(5) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t58 * (pkin(3) + t47) + (-m(3) - m(4)) * t47 - t60 * t45 + t59 * t43) * g(3) + (-m(3) * t40 - m(4) * t55 - t44 * mrSges(2,1) - t46 * mrSges(2,2) - mrSges(1,2) + t58 * (t38 * pkin(7) + t55) + ((m(4) - t58) * qJ(3) + t56) * t39 + t57 * t38) * g(2) + (-m(3) * t41 - m(4) * t52 - t46 * mrSges(2,1) + t44 * mrSges(2,2) - mrSges(1,1) + t58 * (t39 * pkin(7) + t52) + t57 * t39 - t56 * t38) * g(1);
U = t1;
