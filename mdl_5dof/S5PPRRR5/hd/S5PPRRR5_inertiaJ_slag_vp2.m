% Calculate joint inertia matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:39
% EndTime: 2019-12-31 17:35:39
% DurationCPUTime: 0.16s
% Computational Cost: add. (126->54), mult. (287->77), div. (0->0), fcn. (197->6), ass. (0->25)
t23 = cos(qJ(5));
t19 = t23 ^ 2;
t20 = sin(qJ(5));
t35 = t20 ^ 2 + t19;
t39 = mrSges(6,3) * t35;
t21 = sin(qJ(4));
t22 = sin(qJ(3));
t24 = cos(qJ(4));
t25 = cos(qJ(3));
t8 = t21 * t22 - t24 * t25;
t38 = t8 ^ 2;
t36 = Ifges(6,5) * t20 + Ifges(6,6) * t23;
t34 = Ifges(6,2) * t19 + Ifges(5,3) + (Ifges(6,1) * t20 + 0.2e1 * Ifges(6,4) * t23) * t20;
t33 = t35 * pkin(7);
t14 = t21 * pkin(3) + pkin(7);
t32 = t35 * t14;
t31 = -mrSges(6,1) * t20 - mrSges(6,2) * t23;
t30 = 0.2e1 * t39;
t10 = t21 * t25 + t24 * t22;
t11 = -t23 * mrSges(6,1) + t20 * mrSges(6,2);
t29 = (-mrSges(5,1) + t11) * t8 + (-mrSges(5,2) + t39) * t10;
t28 = (t24 * mrSges(5,1) - t21 * mrSges(5,2)) * pkin(3);
t15 = -t24 * pkin(3) - pkin(4);
t5 = t10 ^ 2;
t1 = [m(6) * t35 + m(2) + m(3) + m(4) + m(5); 0; m(3) + m(4) * (t22 ^ 2 + t25 ^ 2) + m(5) * (t5 + t38) + m(6) * (t35 * t5 + t38); 0; t25 * mrSges(4,1) - t22 * mrSges(4,2) + m(6) * (t10 * t32 + t15 * t8) + m(5) * (t10 * t21 - t24 * t8) * pkin(3) + t29; 0.2e1 * t15 * t11 + Ifges(4,3) + 0.2e1 * t28 + t14 * t30 + m(6) * (t35 * t14 ^ 2 + t15 ^ 2) + m(5) * (t21 ^ 2 + t24 ^ 2) * pkin(3) ^ 2 + t34; 0; m(6) * (-pkin(4) * t8 + t10 * t33) + t29; m(6) * (-pkin(4) * t15 + pkin(7) * t32) + (-pkin(4) + t15) * t11 + t28 + (t32 + t33) * mrSges(6,3) + t34; -0.2e1 * pkin(4) * t11 + m(6) * (t35 * pkin(7) ^ 2 + pkin(4) ^ 2) + pkin(7) * t30 + t34; t11; t31 * t10; t31 * t14 + t36; t31 * pkin(7) + t36; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
