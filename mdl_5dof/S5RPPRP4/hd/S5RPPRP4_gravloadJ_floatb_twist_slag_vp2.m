% Calculate Gravitation load on the joints for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:03
% EndTime: 2019-12-31 17:52:04
% DurationCPUTime: 0.23s
% Computational Cost: add. (117->38), mult. (211->43), div. (0->0), fcn. (207->6), ass. (0->20)
t12 = sin(qJ(4));
t13 = cos(qJ(4));
t19 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t23 = -mrSges(5,2) - mrSges(6,2);
t36 = t23 * t12 + t19 * t13;
t35 = m(5) + m(6);
t32 = mrSges(2,1) + mrSges(3,1);
t31 = mrSges(2,2) - mrSges(3,3);
t30 = m(4) + t35;
t29 = t35 * pkin(3) + mrSges(4,1) + t36;
t27 = -m(5) * pkin(6) + m(6) * (-qJ(5) - pkin(6)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t24 = sin(qJ(1));
t25 = cos(qJ(1));
t26 = t25 * pkin(1) + t24 * qJ(2);
t22 = cos(pkin(7));
t21 = sin(pkin(7));
t18 = -t24 * pkin(1) + t25 * qJ(2);
t2 = t25 * t21 - t24 * t22;
t1 = -t24 * t21 - t25 * t22;
t3 = [(-m(3) * t26 - t32 * t25 + t31 * t24 - t30 * (t25 * pkin(2) + t26) + t27 * t2 + t29 * t1) * g(2) + (-m(3) * t18 + t31 * t25 + t32 * t24 - t29 * t2 - t30 * (-t24 * pkin(2) + t18) + t27 * t1) * g(1), (-t24 * g(1) + t25 * g(2)) * (m(3) + t30), t30 * g(3), t36 * g(3) + (g(1) * t1 + g(2) * t2) * (-t19 * t12 + t23 * t13), (-g(1) * t2 + g(2) * t1) * m(6)];
taug = t3(:);
