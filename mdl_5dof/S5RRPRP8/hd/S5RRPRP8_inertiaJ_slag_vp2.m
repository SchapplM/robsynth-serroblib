% Calculate joint inertia matrix for
% S5RRPRP8
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:22
% DurationCPUTime: 0.38s
% Computational Cost: add. (355->127), mult. (625->151), div. (0->0), fcn. (502->4), ass. (0->40)
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t54 = t31 ^ 2 + t33 ^ 2;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t34 = -pkin(2) - pkin(3);
t17 = t32 * qJ(3) + t30 * t34;
t44 = mrSges(5,2) + mrSges(6,2);
t53 = t44 * t17;
t52 = t44 * t30;
t51 = mrSges(5,3) + mrSges(6,3);
t50 = -m(4) * pkin(2) - mrSges(4,1);
t49 = -2 * mrSges(6,1);
t18 = -t33 * pkin(2) - t31 * qJ(3) - pkin(1);
t48 = -0.2e1 * t18;
t47 = (m(6) * pkin(4));
t46 = pkin(6) - pkin(7);
t16 = -t30 * qJ(3) + t32 * t34;
t45 = t16 * mrSges(5,1);
t43 = Ifges(5,5) + Ifges(6,5);
t42 = Ifges(5,6) + Ifges(6,6);
t41 = Ifges(5,3) + Ifges(6,3);
t19 = t46 * t31;
t20 = t46 * t33;
t5 = t30 * t19 + t32 * t20;
t40 = t54 * pkin(6) ^ 2;
t39 = mrSges(6,1) + t47;
t4 = t32 * t19 - t30 * t20;
t6 = t33 * pkin(3) - t18;
t9 = -t31 * t30 - t33 * t32;
t2 = t9 * qJ(5) + t5;
t36 = t4 * mrSges(5,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2);
t15 = t17 ^ 2;
t14 = -pkin(4) + t16;
t10 = -t33 * t30 + t31 * t32;
t8 = t10 * mrSges(6,2);
t7 = t30 * t17;
t3 = -t9 * pkin(4) + t6;
t1 = -t10 * qJ(5) + t4;
t11 = [0.2e1 * t3 * t8 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t48 + (Ifges(4,3) + Ifges(3,2)) * t33) * t33 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t48 + (Ifges(3,1) + Ifges(4,1)) * t31 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t33) * t31 + (-0.2e1 * t6 * mrSges(5,1) + t3 * t49 + 0.2e1 * t5 * mrSges(5,3) + 0.2e1 * t2 * mrSges(6,3) + (Ifges(6,2) + Ifges(5,2)) * t9) * t9 + (0.2e1 * t6 * mrSges(5,2) - 0.2e1 * t4 * mrSges(5,3) - 0.2e1 * t1 * mrSges(6,3) + (Ifges(6,1) + Ifges(5,1)) * t10 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t9) * t10 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(3) * (pkin(1) ^ 2 + t40) + m(4) * (t18 ^ 2 + t40) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(6) * t54; -t1 * mrSges(6,1) + m(6) * (t14 * t1 + t17 * t2) + m(5) * (t16 * t4 + t17 * t5) + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t33 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t31 + (t51 * t17 - t42) * t9 + (-t16 * mrSges(5,3) - t14 * mrSges(6,3) - t43) * t10 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t33 + (-mrSges(3,1) + t50) * t31) * pkin(6) - t36; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t45 + t14 * t49 + 0.2e1 * qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3) + 0.2e1 * t53 + m(5) * (t16 ^ 2 + t15) + m(6) * (t14 ^ 2 + t15) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t41; (m(4) * pkin(6) + mrSges(4,2)) * t31 + m(6) * (t32 * t1 + t30 * t2) + m(5) * (t30 * t5 + t32 * t4) + t51 * (-t32 * t10 + t30 * t9); (-mrSges(5,1) - mrSges(6,1)) * t32 + t52 + m(5) * (t32 * t16 + t7) + m(6) * (t32 * t14 + t7) + t50; m(4) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t30 ^ 2 + t32 ^ 2); t42 * t9 + t39 * t1 + (-mrSges(6,3) * pkin(4) + t43) * t10 + t36; t14 * t47 + t45 - t53 + (-pkin(4) + t14) * mrSges(6,1) - t41; -t52 + (mrSges(5,1) + t39) * t32; (2 * mrSges(6,1) + t47) * pkin(4) + t41; m(6) * t3 - t9 * mrSges(6,1) + t8; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;
