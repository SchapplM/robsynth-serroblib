% Return the minimum parameter vector for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% MPV [17x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t15 = m(6) * pkin(6) + mrSges(6,3);
t6 = sin(pkin(8));
t16 = (mrSges(5,3) + t15) * t6;
t10 = cos(pkin(7));
t7 = sin(pkin(7));
t14 = mrSges(3,1) * t10 - mrSges(3,2) * t7;
t12 = pkin(4) ^ 2;
t5 = sin(pkin(9));
t8 = cos(pkin(9));
t13 = (2 * mrSges(6,3) * pkin(6)) + Ifges(4,1) + Ifges(5,2) + Ifges(6,2) + (-0.2e1 * Ifges(5,4) * t5 + (-m(6) * t12 + Ifges(5,1) - Ifges(5,2)) * t8) * t8;
t11 = pkin(6) ^ 2;
t9 = cos(pkin(8));
t4 = 0.1e1 / t9;
t1 = [(t11 + t12) * m(6) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t6 * ((t15 * pkin(4) - Ifges(5,5)) * t8 + Ifges(5,6) * t5 + Ifges(4,4)) + (-m(6) * t11 + Ifges(4,2) + Ifges(5,3) - t13) * t9) * t9 + 0.2e1 * t14 * pkin(1) + t13; mrSges(2,1) + t14; mrSges(3,1) * t7 + mrSges(3,2) * t10 + mrSges(2,2); m(3); (mrSges(4,1) * t9 - mrSges(4,2) * t6) * t4; mrSges(4,3); m(4); ((m(6) * pkin(4) + mrSges(5,1)) * t9 + t8 * t16) * t4; (mrSges(5,2) * t9 - t5 * t16) * t4; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
