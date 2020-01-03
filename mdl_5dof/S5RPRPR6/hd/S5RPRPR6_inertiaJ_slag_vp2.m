% Calculate joint inertia matrix for
% S5RPRPR6
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (201->64), mult. (336->76), div. (0->0), fcn. (192->6), ass. (0->35)
t24 = cos(qJ(5));
t47 = t24 ^ 2;
t22 = sin(qJ(5));
t36 = t22 ^ 2 + t47;
t35 = t36 * mrSges(6,3);
t10 = -t22 * mrSges(6,1) - t24 * mrSges(6,2);
t46 = mrSges(5,3) - t10;
t21 = cos(pkin(8));
t13 = t21 * pkin(1) + pkin(2);
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t20 = sin(pkin(8));
t43 = pkin(1) * t20;
t7 = t23 * t13 + t25 * t43;
t3 = qJ(4) + t7;
t45 = t3 ^ 2;
t9 = m(6) * t36;
t44 = m(5) + t9;
t6 = t25 * t13 - t23 * t43;
t42 = t6 * mrSges(4,1);
t41 = t7 * mrSges(4,2);
t39 = qJ(4) * t3;
t38 = t24 * mrSges(6,1);
t26 = -pkin(3) - pkin(7);
t34 = t36 * t26;
t33 = 0.2e1 * t46;
t5 = -pkin(3) - t6;
t32 = mrSges(5,2) - t35;
t31 = -t22 * mrSges(6,2) + t38;
t30 = -0.2e1 * t35;
t29 = Ifges(6,1) * t47 + Ifges(5,1) + Ifges(4,3) + (-0.2e1 * Ifges(6,4) * t24 + Ifges(6,2) * t22) * t22;
t27 = qJ(4) ^ 2;
t16 = Ifges(6,5) * t24;
t2 = -pkin(7) + t5;
t1 = [0.2e1 * t42 - 0.2e1 * t41 + 0.2e1 * t5 * mrSges(5,2) + Ifges(2,3) + Ifges(3,3) + t3 * t33 + t2 * t30 + m(5) * (t5 ^ 2 + t45) + m(6) * (t36 * t2 ^ 2 + t45) + m(4) * (t6 ^ 2 + t7 ^ 2) + t29 + (0.2e1 * t21 * mrSges(3,1) - 0.2e1 * t20 * mrSges(3,2) + m(3) * (t20 ^ 2 + t21 ^ 2) * pkin(1)) * pkin(1); 0; m(3) + m(4) + t44; t42 - t41 + (-pkin(3) + t5) * mrSges(5,2) + m(5) * (-pkin(3) * t5 + t39) + m(6) * (t2 * t34 + t39) + t29 + (-t2 - t26) * t35 + t46 * (qJ(4) + t3); 0; -0.2e1 * pkin(3) * mrSges(5,2) + qJ(4) * t33 + t26 * t30 + m(6) * (t36 * t26 ^ 2 + t27) + m(5) * (pkin(3) ^ 2 + t27) + t29; m(5) * t5 + t2 * t9 + t32; 0; -m(5) * pkin(3) + m(6) * t34 + t32; t44; -Ifges(6,6) * t22 + t31 * t2 + t16; t10; t26 * t38 + t16 + (-mrSges(6,2) * t26 - Ifges(6,6)) * t22; t31; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
