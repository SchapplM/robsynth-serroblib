% Calculate joint inertia matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:41
% DurationCPUTime: 0.27s
% Computational Cost: add. (322->93), mult. (452->118), div. (0->0), fcn. (272->6), ass. (0->34)
t30 = cos(qJ(5));
t46 = t30 ^ 2;
t28 = sin(qJ(5));
t40 = t28 ^ 2 + t46;
t39 = mrSges(6,3) * t40;
t45 = -2 * mrSges(5,1);
t44 = 2 * mrSges(5,2);
t43 = 2 * mrSges(4,3);
t13 = t30 * mrSges(6,1) - t28 * mrSges(6,2);
t42 = 0.2e1 * t13;
t31 = cos(qJ(2));
t19 = -t31 * pkin(1) - pkin(2);
t16 = -pkin(3) + t19;
t29 = sin(qJ(2));
t17 = t29 * pkin(1) + qJ(3);
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t5 = t26 * t16 + t27 * t17;
t32 = -pkin(2) - pkin(3);
t10 = t27 * qJ(3) + t26 * t32;
t38 = t26 * t40;
t4 = t27 * t16 - t26 * t17;
t37 = -0.2e1 * t39;
t9 = -t26 * qJ(3) + t27 * t32;
t36 = (t31 * mrSges(3,1) - t29 * mrSges(3,2)) * pkin(1);
t35 = t46 * Ifges(6,2) + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + (Ifges(6,1) * t28 + 0.2e1 * Ifges(6,4) * t30) * t28;
t34 = -mrSges(4,1) + (-mrSges(5,1) - t13) * t27 + (mrSges(5,2) - t39) * t26;
t23 = t27 ^ 2;
t22 = t26 ^ 2;
t8 = -pkin(7) + t10;
t7 = pkin(4) - t9;
t3 = -pkin(7) + t5;
t2 = pkin(4) - t4;
t1 = [-0.2e1 * t19 * mrSges(4,1) + t4 * t45 + t5 * t44 + t17 * t43 + t2 * t42 + Ifges(2,3) + 0.2e1 * t36 + t3 * t37 + m(6) * (t40 * t3 ^ 2 + t2 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(4) * (t17 ^ 2 + t19 ^ 2) + m(3) * (t29 ^ 2 + t31 ^ 2) * pkin(1) ^ 2 + t35; (t2 + t7) * t13 + t36 + (t17 + qJ(3)) * mrSges(4,3) + (t5 + t10) * mrSges(5,2) + (-t4 - t9) * mrSges(5,1) + (-t19 + pkin(2)) * mrSges(4,1) + m(6) * (t40 * t8 * t3 + t7 * t2) + m(5) * (t10 * t5 + t9 * t4) + m(4) * (-pkin(2) * t19 + qJ(3) * t17) + t35 + (-t3 - t8) * t39; 0.2e1 * pkin(2) * mrSges(4,1) + t9 * t45 + t10 * t44 + qJ(3) * t43 + t7 * t42 + t8 * t37 + m(6) * (t40 * t8 ^ 2 + t7 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + m(5) * (t10 ^ 2 + t9 ^ 2) + t35; m(6) * (-t27 * t2 + t3 * t38) + m(5) * (t26 * t5 + t27 * t4) + m(4) * t19 + t34; -m(4) * pkin(2) + m(6) * (-t27 * t7 + t8 * t38) + m(5) * (t26 * t10 + t27 * t9) + t34; m(4) + m(5) * (t22 + t23) + m(6) * (t40 * t22 + t23); 0; 0; 0; m(6) * t40 + m(5); (-mrSges(6,2) * t3 - Ifges(6,6)) * t30 + (-mrSges(6,1) * t3 - Ifges(6,5)) * t28; (-mrSges(6,2) * t8 - Ifges(6,6)) * t30 + (-mrSges(6,1) * t8 - Ifges(6,5)) * t28; (-mrSges(6,1) * t28 - mrSges(6,2) * t30) * t26; t13; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
