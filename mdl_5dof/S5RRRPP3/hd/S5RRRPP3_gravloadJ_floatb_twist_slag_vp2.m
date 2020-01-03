% Calculate Gravitation load on the joints for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:09
% DurationCPUTime: 0.41s
% Computational Cost: add. (261->62), mult. (274->62), div. (0->0), fcn. (218->6), ass. (0->32)
t83 = mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t82 = -mrSges(4,1) + mrSges(5,2);
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t80 = t83 * t30 - t82 * t32;
t29 = qJ(1) + qJ(2);
t24 = sin(t29);
t25 = cos(t29);
t68 = g(1) * t25 + g(2) * t24;
t78 = -t32 * mrSges(6,3) - t80;
t77 = -mrSges(6,1) - mrSges(5,1) - mrSges(4,3) + mrSges(3,2);
t26 = t30 * qJ(4);
t27 = t32 * pkin(3);
t54 = t27 + t26;
t70 = m(5) + m(6);
t44 = m(6) * (-pkin(3) - qJ(5)) - mrSges(6,3);
t47 = -pkin(2) - t26;
t66 = t77 * t25 + (-m(6) * t47 - t44 * t32 - m(5) * (t47 - t27) + mrSges(3,1) + t80) * t24;
t65 = (-mrSges(3,1) + t78) * t25 + t77 * t24;
t31 = sin(qJ(1));
t62 = t31 * pkin(1);
t33 = cos(qJ(1));
t28 = t33 * pkin(1);
t55 = t25 * pkin(2) + t24 * pkin(7);
t52 = t32 * qJ(5);
t21 = t25 * pkin(7);
t50 = t21 - t62;
t49 = -t24 * pkin(2) + t21;
t46 = t54 * t25 + t55;
t39 = t24 * pkin(4) + t25 * t52 + t46;
t22 = t25 * pkin(4);
t1 = [(-t33 * mrSges(2,1) + t31 * mrSges(2,2) - m(3) * t28 - m(4) * (t28 + t55) - m(5) * (t28 + t46) - m(6) * (t28 + t39) + t65) * g(2) + (t31 * mrSges(2,1) + t33 * mrSges(2,2) + m(3) * t62 - m(4) * (t49 - t62) - m(5) * t50 - m(6) * (t22 + t50) + t66) * g(1), (-m(4) * t55 - m(5) * t46 - m(6) * t39 + t65) * g(2) + (-m(4) * t49 - m(5) * t21 - m(6) * (t21 + t22) + t66) * g(1), (-m(5) * t54 - m(6) * (t52 + t54) + t78) * g(3) + ((m(5) * pkin(3) - t44 - t82) * t30 + (-qJ(4) * t70 - t83) * t32) * t68, (g(3) * t32 - t30 * t68) * t70, (-g(3) * t30 - t32 * t68) * m(6)];
taug = t1(:);
