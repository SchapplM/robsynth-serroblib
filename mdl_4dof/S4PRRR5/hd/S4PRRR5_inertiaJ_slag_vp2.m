% Calculate joint inertia matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:36
% DurationCPUTime: 0.17s
% Computational Cost: add. (120->50), mult. (279->76), div. (0->0), fcn. (193->6), ass. (0->25)
t23 = cos(qJ(4));
t19 = t23 ^ 2;
t20 = sin(qJ(4));
t35 = t20 ^ 2 + t19;
t39 = mrSges(5,3) * t35;
t21 = sin(qJ(3));
t22 = sin(qJ(2));
t24 = cos(qJ(3));
t25 = cos(qJ(2));
t8 = t21 * t22 - t24 * t25;
t38 = t8 ^ 2;
t36 = Ifges(5,5) * t20 + Ifges(5,6) * t23;
t34 = Ifges(5,2) * t19 + Ifges(4,3) + (Ifges(5,1) * t20 + 0.2e1 * Ifges(5,4) * t23) * t20;
t33 = t35 * pkin(6);
t14 = t21 * pkin(2) + pkin(6);
t32 = t35 * t14;
t31 = -mrSges(5,1) * t20 - mrSges(5,2) * t23;
t30 = 0.2e1 * t39;
t10 = t21 * t25 + t24 * t22;
t11 = -t23 * mrSges(5,1) + t20 * mrSges(5,2);
t29 = (-mrSges(4,1) + t11) * t8 + (-mrSges(4,2) + t39) * t10;
t28 = (t24 * mrSges(4,1) - t21 * mrSges(4,2)) * pkin(2);
t15 = -t24 * pkin(2) - pkin(3);
t5 = t10 ^ 2;
t1 = [m(2) + m(3) * (t22 ^ 2 + t25 ^ 2) + m(4) * (t5 + t38) + m(5) * (t35 * t5 + t38); t25 * mrSges(3,1) - t22 * mrSges(3,2) + m(5) * (t10 * t32 + t15 * t8) + m(4) * (t10 * t21 - t24 * t8) * pkin(2) + t29; 0.2e1 * t15 * t11 + Ifges(3,3) + 0.2e1 * t28 + t14 * t30 + m(5) * (t35 * t14 ^ 2 + t15 ^ 2) + m(4) * (t21 ^ 2 + t24 ^ 2) * pkin(2) ^ 2 + t34; m(5) * (-pkin(3) * t8 + t10 * t33) + t29; m(5) * (-pkin(3) * t15 + pkin(6) * t32) + (t15 - pkin(3)) * t11 + t28 + (t32 + t33) * mrSges(5,3) + t34; -0.2e1 * pkin(3) * t11 + m(5) * (t35 * pkin(6) ^ 2 + pkin(3) ^ 2) + pkin(6) * t30 + t34; t31 * t10; t31 * t14 + t36; t31 * pkin(6) + t36; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
