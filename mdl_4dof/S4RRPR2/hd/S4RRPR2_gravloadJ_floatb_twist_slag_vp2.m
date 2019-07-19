% Calculate Gravitation load on the joints for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:30
% EndTime: 2019-07-18 18:16:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (134->35), mult. (116->38), div. (0->0), fcn. (98->6), ass. (0->20)
t36 = -mrSges(3,1) - mrSges(4,1);
t35 = mrSges(3,2) - mrSges(4,3);
t34 = m(4) + m(5);
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t19 = cos(t21);
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t5 = -t18 * t29 - t19 * t30;
t6 = -t18 * t30 + t19 * t29;
t33 = -t6 * mrSges(5,1) + t5 * mrSges(5,2);
t32 = t5 * mrSges(5,1) + t6 * mrSges(5,2);
t28 = t19 * pkin(2) + t18 * qJ(3);
t23 = cos(qJ(1));
t27 = t23 * pkin(1) + t28;
t25 = t35 * t18 + t36 * t19 + t32;
t24 = (m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(3)) - t36) * t18 + t33 + (-qJ(3) * t34 + t35) * t19;
t22 = sin(qJ(1));
t16 = t19 * pkin(3);
t1 = [(t22 * mrSges(2,2) - m(4) * t27 - m(5) * (t16 + t27) + (-m(3) * pkin(1) - mrSges(2,1)) * t23 + t25) * g(2) + (mrSges(2,2) * t23 + (mrSges(2,1) + (m(3) + t34) * pkin(1)) * t22 + t24) * g(1), (-m(4) * t28 - m(5) * (t16 + t28) + t25) * g(2) + t24 * g(1), t34 * (-g(1) * t18 + g(2) * t19), -g(1) * t33 - g(2) * t32];
taug  = t1(:);
