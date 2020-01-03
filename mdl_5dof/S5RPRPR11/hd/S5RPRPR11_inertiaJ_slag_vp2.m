% Calculate joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:12
% DurationCPUTime: 0.28s
% Computational Cost: add. (436->111), mult. (808->145), div. (0->0), fcn. (797->6), ass. (0->40)
t29 = cos(pkin(8));
t27 = t29 ^ 2;
t24 = -t29 * pkin(2) - pkin(1);
t48 = 0.2e1 * t24;
t33 = -pkin(3) - pkin(4);
t47 = cos(qJ(3));
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t20 = -qJ(4) * t30 + t32 * t33;
t46 = t20 * mrSges(6,1);
t21 = qJ(4) * t32 + t30 * t33;
t45 = t21 * mrSges(6,2);
t44 = mrSges(5,2) + mrSges(4,3);
t43 = pkin(6) + qJ(2);
t22 = t43 * t29;
t31 = sin(qJ(3));
t28 = sin(pkin(8));
t40 = t43 * t28;
t12 = t47 * t22 - t31 * t40;
t42 = t28 ^ 2 + t27;
t10 = t22 * t31 + t47 * t40;
t41 = t10 ^ 2 + t12 ^ 2;
t39 = -t29 * mrSges(3,1) + t28 * mrSges(3,2);
t18 = t31 * t28 - t29 * t47;
t19 = t28 * t47 + t31 * t29;
t7 = t18 * t32 - t19 * t30;
t8 = t18 * t30 + t19 * t32;
t38 = -t7 * mrSges(6,1) + t8 * mrSges(6,2);
t37 = t32 * mrSges(6,1) - t30 * mrSges(6,2);
t36 = t19 * qJ(4) - t24;
t4 = -pkin(7) * t19 + t10;
t5 = pkin(7) * t18 + t12;
t1 = -t30 * t5 + t32 * t4;
t2 = t30 * t4 + t32 * t5;
t35 = t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t8 + Ifges(6,6) * t7;
t14 = t19 * mrSges(4,2);
t13 = t18 * mrSges(5,1);
t6 = pkin(3) * t18 - t36;
t3 = t18 * t33 + t36;
t9 = [Ifges(3,2) * t27 - 0.2e1 * pkin(1) * t39 + 0.2e1 * t6 * t13 + t14 * t48 + 0.2e1 * t3 * t38 + t8 * (Ifges(6,1) * t8 + Ifges(6,4) * t7) + t7 * (Ifges(6,4) * t8 + Ifges(6,2) * t7) + Ifges(2,3) + (Ifges(3,1) * t28 + 0.2e1 * Ifges(3,4) * t29) * t28 + 0.2e1 * (-t1 * t8 + t2 * t7) * mrSges(6,3) + 0.2e1 * t42 * qJ(2) * mrSges(3,3) + (-0.2e1 * t6 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t19 + 0.2e1 * t44 * t10) * t19 + (mrSges(4,1) * t48 + (Ifges(5,3) + Ifges(4,2)) * t18 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t19 - 0.2e1 * t44 * t12) * t18 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t6 ^ 2 + t41) + m(4) * (t24 ^ 2 + t41) + m(3) * (qJ(2) ^ 2 * t42 + pkin(1) ^ 2); -m(3) * pkin(1) + m(4) * t24 + m(5) * t6 - m(6) * t3 + t18 * mrSges(4,1) - t19 * mrSges(5,3) + t13 + t14 - t38 + t39; m(3) + m(4) + m(5) + m(6); (-mrSges(4,2) + mrSges(5,3)) * t12 + (-mrSges(4,1) - mrSges(5,1)) * t10 + (-t20 * t8 + t21 * t7) * mrSges(6,3) + m(6) * (t1 * t20 + t2 * t21) + m(5) * (-pkin(3) * t10 + qJ(4) * t12) + (-mrSges(5,2) * pkin(3) + Ifges(5,4) + Ifges(4,5)) * t19 + (-mrSges(5,2) * qJ(4) - Ifges(4,6) + Ifges(5,6)) * t18 - t35; 0; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t46 + 0.2e1 * t45 + 0.2e1 * qJ(4) * mrSges(5,3) + Ifges(5,2) + Ifges(4,3) + Ifges(6,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2); t19 * mrSges(5,2) + (t30 * t7 - t32 * t8) * mrSges(6,3) + m(6) * (t1 * t32 + t2 * t30) + m(5) * t10; 0; -mrSges(5,1) - m(5) * pkin(3) + m(6) * (t20 * t32 + t21 * t30) - t37; m(5) + m(6) * (t30 ^ 2 + t32 ^ 2); t35; 0; -Ifges(6,3) - t45 + t46; t37; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
