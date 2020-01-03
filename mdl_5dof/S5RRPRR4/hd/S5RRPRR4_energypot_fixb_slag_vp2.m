% Calculate potential energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:46
% EndTime: 2020-01-03 12:01:46
% DurationCPUTime: 0.24s
% Computational Cost: add. (136->49), mult. (92->40), div. (0->0), fcn. (61->10), ass. (0->20)
t57 = -m(4) - m(5) - m(6);
t42 = qJ(4) + qJ(5);
t36 = sin(t42);
t38 = cos(t42);
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t56 = mrSges(4,1) + m(5) * pkin(3) + t46 * mrSges(5,1) - t44 * mrSges(5,2) + m(6) * (t46 * pkin(4) + pkin(3)) + t38 * mrSges(6,1) - t36 * mrSges(6,2);
t55 = m(5) * pkin(7) - m(6) * (-pkin(8) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t54 = pkin(6) + pkin(5);
t45 = sin(qJ(1));
t40 = t45 * pkin(1);
t47 = cos(qJ(1));
t53 = t47 * pkin(1);
t43 = qJ(1) + qJ(2);
t39 = cos(t43);
t37 = sin(t43);
t35 = pkin(9) + t43;
t32 = cos(t35);
t31 = sin(t35);
t1 = (m(3) * t53 + t47 * mrSges(2,1) + t39 * mrSges(3,1) - t45 * mrSges(2,2) - t37 * mrSges(3,2) - mrSges(1,3) + t57 * (-pkin(2) * t39 - t53) + t56 * t32 + t55 * t31) * g(3) + (-m(3) * t40 - t45 * mrSges(2,1) - t37 * mrSges(3,1) - t47 * mrSges(2,2) - t39 * mrSges(3,2) - mrSges(1,2) + t57 * (pkin(2) * t37 + t40) + t55 * t32 - t56 * t31) * g(2) + (-m(2) * pkin(5) - m(3) * t54 - t36 * mrSges(6,1) - t46 * mrSges(5,2) - t38 * mrSges(6,2) - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t44 + t57 * (qJ(3) + t54)) * g(1);
U = t1;
