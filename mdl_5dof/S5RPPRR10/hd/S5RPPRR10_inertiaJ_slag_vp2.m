% Calculate joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:01
% DurationCPUTime: 0.36s
% Computational Cost: add. (372->102), mult. (727->146), div. (0->0), fcn. (710->6), ass. (0->35)
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t48 = t30 ^ 2 + t31 ^ 2;
t43 = t30 * qJ(3) + pkin(1);
t22 = -t31 * pkin(2) - t43;
t47 = -0.2e1 * t22;
t46 = -pkin(6) + qJ(2);
t23 = t46 * t30;
t24 = t46 * t31;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t12 = t33 * t23 + t35 * t24;
t45 = t48 * qJ(2) ^ 2;
t32 = sin(qJ(5));
t34 = cos(qJ(5));
t20 = -t32 * t33 + t34 * t35;
t21 = t32 * t35 + t34 * t33;
t44 = t20 * mrSges(6,1) - t21 * mrSges(6,2);
t11 = t35 * t23 - t33 * t24;
t18 = -t30 * t33 - t31 * t35;
t19 = t30 * t35 - t31 * t33;
t8 = t34 * t18 - t32 * t19;
t9 = t32 * t18 + t34 * t19;
t41 = -t8 * mrSges(6,1) + t9 * mrSges(6,2);
t4 = -t19 * pkin(7) + t11;
t5 = t18 * pkin(7) + t12;
t2 = -t32 * t5 + t34 * t4;
t3 = t32 * t4 + t34 * t5;
t40 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t9 + Ifges(6,6) * t8;
t39 = -t18 * mrSges(5,1) + t19 * mrSges(5,2);
t38 = (t34 * mrSges(6,1) - t32 * mrSges(6,2)) * pkin(4);
t13 = (pkin(2) + pkin(3)) * t31 + t43;
t27 = t30 * mrSges(3,2);
t10 = -t18 * pkin(4) + t13;
t1 = [Ifges(2,3) + Ifges(6,1) * t9 ^ 2 + 0.2e1 * t10 * t41 + Ifges(5,1) * t19 ^ 2 + 0.2e1 * t13 * t39 - 0.2e1 * pkin(1) * t27 + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t47 + (Ifges(4,3) + Ifges(3,2)) * t31) * t31 + m(6) * (t10 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2 + t13 ^ 2) + m(4) * (t22 ^ 2 + t45) + m(3) * (pkin(1) ^ 2 + t45) + (0.2e1 * Ifges(6,4) * t9 + Ifges(6,2) * t8) * t8 + (0.2e1 * Ifges(5,4) * t19 + Ifges(5,2) * t18) * t18 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * qJ(2) * t48 + 0.2e1 * (-t2 * t9 + t3 * t8) * mrSges(6,3) + 0.2e1 * (-t11 * t19 + t12 * t18) * mrSges(5,3) + (mrSges(4,3) * t47 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t31 + (Ifges(4,1) + Ifges(3,1)) * t30) * t30; -m(3) * pkin(1) - t30 * mrSges(4,3) + t27 + (-mrSges(4,1) - mrSges(3,1)) * t31 + m(4) * t22 - m(5) * t13 - m(6) * t10 - t39 - t41; m(3) + m(4) + m(5) + m(6); (m(4) * qJ(2) + mrSges(4,2)) * t30 + (-t20 * t9 + t21 * t8) * mrSges(6,3) + (t33 * t18 - t35 * t19) * mrSges(5,3) + m(6) * (t20 * t2 + t21 * t3) + m(5) * (t35 * t11 + t33 * t12); 0; m(4) + m(5) * (t33 ^ 2 + t35 ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2); t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,5) * t19 + Ifges(5,6) * t18 + (m(6) * (t2 * t34 + t3 * t32) + (t32 * t8 - t34 * t9) * mrSges(6,3)) * pkin(4) + t40; 0; t35 * mrSges(5,1) - t33 * mrSges(5,2) + m(6) * (t20 * t34 + t21 * t32) * pkin(4) + t44; Ifges(5,3) + Ifges(6,3) + m(6) * (t32 ^ 2 + t34 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t38; t40; 0; t44; Ifges(6,3) + t38; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
