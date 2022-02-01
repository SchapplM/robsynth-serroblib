% Return the minimum parameter vector for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t13 = (m(6) * pkin(8));
t7 = pkin(4) ^ 2 * m(6);
t2 = -Ifges(5,1) + Ifges(5,2) + t7;
t5 = cos(pkin(9));
t3 = t5 ^ 2;
t12 = t2 * t3;
t4 = sin(pkin(9));
t11 = t4 * t5;
t10 = Ifges(5,4) * t11;
t9 = -mrSges(6,3) - t13;
t6 = m(3) + m(4);
t1 = t9 * pkin(4) + Ifges(5,5);
t8 = [t6 * pkin(1) ^ 2 + Ifges(2,3); t6 * pkin(1) + mrSges(2,1); mrSges(2,2); t12 + 0.2e1 * t10 + Ifges(5,1) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)) + (2 * pkin(7) * mrSges(4,3)) + Ifges(4,2) + Ifges(6,2) + Ifges(3,3) + ((2 * mrSges(6,3) + t13) * pkin(8)); pkin(2) * m(4) + mrSges(3,1); -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3); Ifges(4,1) - Ifges(4,2) - 0.4e1 * t10 + t2 - 0.2e1 * t12; 0.2e1 * Ifges(5,4) * t3 - t2 * t11 + Ifges(4,4) - Ifges(5,4); -t4 * Ifges(5,6) + t1 * t5 + Ifges(4,5); t5 * Ifges(5,6) + t1 * t4 + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + t7; mrSges(4,1); mrSges(4,2); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t9; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t8;
