% Calculate joint inertia matrix for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:49
% EndTime: 2019-03-08 18:31:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (107->37), mult. (198->47), div. (0->0), fcn. (138->6), ass. (0->23)
t10 = sin(pkin(7));
t25 = pkin(1) * t10;
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t11 = cos(pkin(7));
t8 = t11 * pkin(1) + pkin(2);
t5 = -t13 * t25 + t15 * t8;
t4 = pkin(3) + t5;
t6 = t13 * t8 + t15 * t25;
t3 = t12 * t4 + t14 * t6;
t24 = t3 * mrSges(5,2);
t23 = t5 * mrSges(4,1);
t22 = t6 * mrSges(4,2);
t21 = t12 * mrSges(5,2);
t20 = Ifges(4,3) + Ifges(5,3);
t19 = pkin(3) * t21;
t2 = -t12 * t6 + t14 * t4;
t1 = t2 * mrSges(5,1);
t18 = Ifges(5,3) + t1 - t24;
t9 = t14 * pkin(3) * mrSges(5,1);
t7 = [0.2e1 * t23 - 0.2e1 * t22 - 0.2e1 * t24 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t1 + m(5) * (t2 ^ 2 + t3 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + t20 + (0.2e1 * t11 * mrSges(3,1) - 0.2e1 * t10 * mrSges(3,2) + m(3) * (t10 ^ 2 + t11 ^ 2) * pkin(1)) * pkin(1); 0; m(3) + m(4) + m(5); t23 - t22 + Ifges(4,3) + t9 + (m(5) * (t12 * t3 + t14 * t2) - t21) * pkin(3) + t18; 0; -0.2e1 * t19 + 0.2e1 * t9 + m(5) * (t12 ^ 2 + t14 ^ 2) * pkin(3) ^ 2 + t20; t18; 0; Ifges(5,3) + t9 - t19; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
