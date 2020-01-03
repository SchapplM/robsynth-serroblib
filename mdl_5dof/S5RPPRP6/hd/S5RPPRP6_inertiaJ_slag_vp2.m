% Calculate joint inertia matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:54:59
% DurationCPUTime: 0.26s
% Computational Cost: add. (261->77), mult. (446->93), div. (0->0), fcn. (372->4), ass. (0->33)
t46 = m(5) + m(6);
t47 = (-mrSges(6,2) - mrSges(5,3));
t24 = sin(pkin(7));
t27 = sin(qJ(4));
t25 = cos(pkin(7));
t39 = cos(qJ(4));
t31 = t39 * t25;
t12 = t27 * t24 - t31;
t38 = t27 * t25;
t13 = t24 * t39 + t38;
t34 = t12 ^ 2 + t13 ^ 2;
t45 = -m(6) * pkin(4) - mrSges(6,1);
t44 = -mrSges(5,1) + t45;
t43 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t22 = t25 ^ 2;
t42 = 2 * mrSges(6,3);
t26 = -pkin(1) - qJ(3);
t40 = -pkin(6) + t26;
t37 = t24 * mrSges(4,1) + t25 * mrSges(4,2);
t36 = t24 ^ 2 + t22;
t16 = t24 * pkin(3) + qJ(2);
t15 = t40 * t24;
t4 = t27 * t15 - t31 * t40;
t6 = t15 * t39 + t38 * t40;
t35 = t4 ^ 2 + t6 ^ 2;
t32 = m(4) * t36;
t30 = t36 * mrSges(4,3);
t29 = 2 * t47;
t28 = qJ(2) ^ 2;
t9 = t12 * mrSges(5,2);
t8 = t13 * mrSges(6,1);
t2 = t13 * pkin(4) + t12 * qJ(5) + t16;
t1 = [Ifges(3,1) + Ifges(2,3) - 0.2e1 * t16 * t9 + 0.2e1 * t2 * t8 - (2 * pkin(1) * mrSges(3,2)) + Ifges(4,1) * t22 + (-0.2e1 * Ifges(4,4) * t25 + Ifges(4,2) * t24) * t24 + (0.2e1 * t16 * mrSges(5,1) + (Ifges(6,3) + Ifges(5,2)) * t13 + t6 * t29) * t13 + m(5) * (t16 ^ 2 + t35) + m(6) * (t2 ^ 2 + t35) + m(4) * (t26 ^ 2 * t36 + t28) + m(3) * ((pkin(1) ^ 2) + t28) - 0.2e1 * t26 * t30 + 0.2e1 * (mrSges(3,3) + t37) * qJ(2) + (t2 * t42 + t4 * t29 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t13 + (Ifges(5,1) + Ifges(6,1)) * t12) * t12; -m(3) * pkin(1) + mrSges(3,2) - t30 + t26 * t32 + t46 * (t12 * t4 + t13 * t6) + t47 * t34; t46 * t34 + m(3) + t32; m(4) * qJ(2) + m(5) * t16 + m(6) * t2 + t13 * mrSges(5,1) + t12 * mrSges(6,3) + t37 + t8 - t9; 0; m(4) + t46; (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t13 + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t12 + t43 * t6 + t44 * t4; t44 * t12 + t43 * t13; 0; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + qJ(5) * t42 + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); m(6) * t4 - t12 * mrSges(6,2); m(6) * t12; 0; t45; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
