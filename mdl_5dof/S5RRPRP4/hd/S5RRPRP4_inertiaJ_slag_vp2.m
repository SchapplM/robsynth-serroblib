% Calculate joint inertia matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:38
% DurationCPUTime: 0.30s
% Computational Cost: add. (270->86), mult. (443->93), div. (0->0), fcn. (201->4), ass. (0->35)
t33 = cos(qJ(4));
t29 = t33 ^ 2;
t31 = sin(qJ(4));
t51 = t31 ^ 2 + t29;
t60 = (mrSges(5,3) + mrSges(6,2)) * t51;
t62 = t31 * mrSges(5,1) + t33 * mrSges(5,2) + mrSges(4,3);
t32 = sin(qJ(2));
t55 = t32 * pkin(1);
t18 = qJ(3) + t55;
t59 = t18 ^ 2;
t12 = t31 * mrSges(6,1) - t33 * mrSges(6,3);
t58 = 0.2e1 * t12;
t34 = cos(qJ(2));
t20 = -t34 * pkin(1) - pkin(2);
t17 = -pkin(7) + t20;
t35 = -pkin(2) - pkin(7);
t57 = t51 * t17 * t35;
t56 = m(6) * t33;
t53 = t51 * t17 ^ 2;
t52 = t51 * t35 ^ 2;
t50 = qJ(3) * t18;
t48 = 0.2e1 * t62;
t9 = t31 * pkin(4) - t33 * qJ(5) + qJ(3);
t46 = pkin(4) * t33 + qJ(5) * t31;
t45 = (t34 * mrSges(3,1) - t32 * mrSges(3,2)) * pkin(1);
t44 = 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t51;
t43 = Ifges(4,1) + Ifges(3,3) + (Ifges(6,1) + Ifges(5,1)) * t29 + ((Ifges(5,2) + Ifges(6,3)) * t31 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t33) * t31;
t41 = mrSges(4,2) - t60;
t40 = 0.2e1 * t60;
t39 = -t46 * mrSges(6,2) + (Ifges(6,4) + Ifges(5,5)) * t33 + (-Ifges(5,6) + Ifges(6,6)) * t31;
t38 = m(6) * t46 + (mrSges(5,1) + mrSges(6,1)) * t33 + (-mrSges(5,2) + mrSges(6,3)) * t31;
t36 = qJ(3) ^ 2;
t23 = t33 * mrSges(6,2);
t2 = t9 + t55;
t1 = [0.2e1 * t20 * mrSges(4,2) + t2 * t58 + Ifges(2,3) + t18 * t48 + 0.2e1 * t45 - t40 * t17 + m(6) * (t2 ^ 2 + t53) + m(5) * (t53 + t59) + m(4) * (t20 ^ 2 + t59) + m(3) * (t32 ^ 2 + t34 ^ 2) * pkin(1) ^ 2 + t43; (t9 + t2) * t12 + t45 + (t20 - pkin(2)) * mrSges(4,2) + m(6) * (t9 * t2 + t57) + m(5) * (t50 + t57) + m(4) * (-pkin(2) * t20 + t50) + t43 + t62 * (t18 + qJ(3)) + (-t17 - t35) * t60; -0.2e1 * pkin(2) * mrSges(4,2) + t9 * t58 + qJ(3) * t48 + m(6) * (t9 ^ 2 + t52) + m(5) * (t36 + t52) + m(4) * (pkin(2) ^ 2 + t36) - t40 * t35 + t43; m(4) * t20 + t17 * t44 + t41; -m(4) * pkin(2) + t35 * t44 + t41; m(4) + t44; t38 * t17 + t39; t38 * t35 + t39; t38; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); -t17 * t56 + t23; -t35 * t56 + t23; -t56; -m(6) * pkin(4) - mrSges(6,1); m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
