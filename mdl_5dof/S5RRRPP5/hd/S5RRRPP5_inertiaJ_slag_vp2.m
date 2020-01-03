% Calculate joint inertia matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:37
% DurationCPUTime: 0.36s
% Computational Cost: add. (406->131), mult. (739->157), div. (0->0), fcn. (603->4), ass. (0->39)
t56 = (mrSges(5,3) + mrSges(6,2));
t55 = 2 * mrSges(5,1);
t54 = -2 * mrSges(6,1);
t35 = cos(qJ(2));
t29 = -t35 * pkin(2) - pkin(1);
t53 = 0.2e1 * t29;
t36 = -pkin(3) - pkin(4);
t52 = -pkin(7) - pkin(6);
t51 = -mrSges(6,1) - mrSges(5,1);
t50 = mrSges(5,2) + mrSges(4,3);
t49 = mrSges(5,2) - mrSges(6,3);
t22 = t52 * t35;
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t33 = sin(qJ(2));
t44 = t52 * t33;
t11 = -t34 * t22 + t32 * t44;
t48 = t33 ^ 2 + t35 ^ 2;
t9 = -t32 * t22 - t34 * t44;
t46 = t11 ^ 2 + t9 ^ 2;
t45 = Ifges(5,2) + Ifges(4,3) + Ifges(6,3);
t28 = -t34 * pkin(2) - pkin(3);
t43 = 2 * t56;
t18 = t32 * t35 + t34 * t33;
t42 = t18 * qJ(4) - t29;
t41 = (t34 * mrSges(4,1) - t32 * mrSges(4,2)) * pkin(2);
t17 = t32 * t33 - t34 * t35;
t3 = -t18 * qJ(5) + t9;
t4 = t17 * qJ(5) + t11;
t40 = -t3 * mrSges(6,1) + t4 * mrSges(6,2) + (-mrSges(4,1) - mrSges(5,1)) * t9 + (Ifges(5,4) + Ifges(4,5)) * t18 + (-Ifges(4,6) + Ifges(5,6)) * t17 + (-mrSges(4,2) + mrSges(5,3)) * t11;
t37 = qJ(4) ^ 2;
t26 = t32 * pkin(2) + qJ(4);
t25 = -pkin(4) + t28;
t24 = t26 ^ 2;
t23 = qJ(4) * t26;
t12 = t18 * mrSges(6,2);
t5 = t17 * pkin(3) - t42;
t1 = t36 * t17 + t42;
t2 = [t33 * (Ifges(3,1) * t33 + Ifges(3,4) * t35) + t35 * (Ifges(3,4) * t33 + Ifges(3,2) * t35) - 0.2e1 * pkin(1) * (-t35 * mrSges(3,1) + t33 * mrSges(3,2)) + 0.2e1 * t1 * t12 + Ifges(2,3) + 0.2e1 * t48 * pkin(6) * mrSges(3,3) + m(6) * (t1 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(5) * (t5 ^ 2 + t46) + m(4) * (t29 ^ 2 + t46) + m(3) * (t48 * pkin(6) ^ 2 + pkin(1) ^ 2) + (mrSges(4,1) * t53 + t5 * t55 + t1 * t54 + 0.2e1 * t4 * mrSges(6,3) + (Ifges(5,3) + Ifges(6,2) + Ifges(4,2)) * t17 - 0.2e1 * t50 * t11) * t17 + (mrSges(4,2) * t53 - 0.2e1 * t5 * mrSges(5,3) - 0.2e1 * t3 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1) + Ifges(4,1)) * t18 + 0.2e1 * t50 * t9 + 0.2e1 * (-Ifges(4,4) + Ifges(6,4) + Ifges(5,5)) * t17) * t18; t40 + (m(4) * (t11 * t32 - t34 * t9) + (-t32 * t17 - t34 * t18) * mrSges(4,3)) * pkin(2) + (-t49 * t26 - Ifges(6,6)) * t17 + (t28 * mrSges(5,2) - t25 * mrSges(6,3) - Ifges(6,5)) * t18 + m(6) * (t25 * t3 + t26 * t4) + m(5) * (t26 * t11 + t28 * t9) + (-t33 * mrSges(3,1) - t35 * mrSges(3,2)) * pkin(6) + Ifges(3,5) * t33 + Ifges(3,6) * t35; -0.2e1 * t28 * mrSges(5,1) + t25 * t54 + Ifges(3,3) + t26 * t43 + 0.2e1 * t41 + m(5) * (t28 ^ 2 + t24) + m(6) * (t25 ^ 2 + t24) + m(4) * (t32 ^ 2 + t34 ^ 2) * pkin(2) ^ 2 + t45; t40 + (-t49 * qJ(4) - Ifges(6,6)) * t17 + m(5) * (-pkin(3) * t9 + qJ(4) * t11) + m(6) * (qJ(4) * t4 + t36 * t3) + (-pkin(3) * mrSges(5,2) - t36 * mrSges(6,3) - Ifges(6,5)) * t18; t41 + (-t36 - t25) * mrSges(6,1) + (-t28 + pkin(3)) * mrSges(5,1) + m(5) * (-pkin(3) * t28 + t23) + m(6) * (t36 * t25 + t23) + t45 + t56 * (t26 + qJ(4)); pkin(3) * t55 + t36 * t54 + qJ(4) * t43 + m(5) * (pkin(3) ^ 2 + t37) + m(6) * (t36 ^ 2 + t37) + t45; m(5) * t9 + m(6) * t3 + t49 * t18; m(5) * t28 + m(6) * t25 + t51; -m(5) * pkin(3) + m(6) * t36 + t51; m(5) + m(6); m(6) * t1 - t17 * mrSges(6,1) + t12; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
