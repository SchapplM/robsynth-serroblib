% Return the minimum parameter vector for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t5 = m(5) + m(6);
t2 = m(4) + t5;
t1 = t2 * pkin(2) + mrSges(3,1);
t3 = sin(pkin(9));
t4 = cos(pkin(9));
t6 = -mrSges(3,2) * t3 + t1 * t4;
t7 = [t2 * pkin(2) ^ 2 + 0.2e1 * pkin(1) * t6 + Ifges(2,3) + Ifges(3,3); mrSges(2,1) + t6; mrSges(3,2) * t4 + t1 * t3 + mrSges(2,2); m(3) + t2; t5 * pkin(3) ^ 2 + Ifges(4,3); t5 * pkin(3) + mrSges(4,1); mrSges(4,2); Ifges(5,3) + Ifges(6,2) + 2 * pkin(8) * mrSges(6,3) + (pkin(4) ^ 2 + pkin(8) ^ 2) * m(6); pkin(4) * m(6) + mrSges(5,1); -pkin(8) * m(6) + mrSges(5,2) - mrSges(6,3); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t7;
