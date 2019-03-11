% Calculate joint inertia matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:26
% EndTime: 2019-03-08 18:28:27
% DurationCPUTime: 0.09s
% Computational Cost: add. (100->39), mult. (139->50), div. (0->0), fcn. (98->4), ass. (0->16)
t11 = sin(qJ(4));
t12 = cos(qJ(4));
t10 = cos(pkin(6));
t13 = -pkin(1) - pkin(2);
t9 = sin(pkin(6));
t6 = -t9 * qJ(2) + t10 * t13;
t5 = -pkin(3) + t6;
t7 = t10 * qJ(2) + t9 * t13;
t1 = -t11 * t7 + t12 * t5;
t16 = t1 * mrSges(5,1);
t2 = t11 * t5 + t12 * t7;
t15 = t2 * mrSges(5,2);
t3 = t12 * t10 - t11 * t9;
t4 = t11 * t10 + t12 * t9;
t14 = t3 * mrSges(5,1) - t4 * mrSges(5,2);
t8 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t6 * mrSges(4,1) - 0.2e1 * t16 + 0.2e1 * t7 * mrSges(4,2) + 0.2e1 * t15 + 0.2e1 * qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + Ifges(5,3) + m(5) * (t1 ^ 2 + t2 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + m(4) * (t6 ^ 2 + t7 ^ 2); -m(3) * pkin(1) - t10 * mrSges(4,1) + t9 * mrSges(4,2) - mrSges(3,1) + m(5) * (t3 * t1 + t4 * t2) + m(4) * (t10 * t6 + t9 * t7) - t14; m(3) + m(4) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2); 0; 0; m(4) + m(5); -Ifges(5,3) - t15 + t16; t14; 0; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1) t8(2) t8(4) t8(7); t8(2) t8(3) t8(5) t8(8); t8(4) t8(5) t8(6) t8(9); t8(7) t8(8) t8(9) t8(10);];
Mq  = res;
