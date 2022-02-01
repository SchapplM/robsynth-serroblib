% Return the minimum parameter vector for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t6 = (pkin(7) * m(6));
t2 = sin(pkin(8));
t4 = cos(pkin(8));
t5 = mrSges(4,1) * t4 - mrSges(4,2) * t2;
t3 = cos(pkin(9));
t1 = sin(pkin(9));
t7 = [pkin(1) ^ 2 * m(3) + Ifges(2,3); pkin(1) * m(3) + mrSges(2,1); mrSges(2,2); Ifges(6,2) + Ifges(3,3) + Ifges(4,3) + Ifges(5,1) + (0.2e1 * t1 * Ifges(5,4) + (m(6) * pkin(4) ^ 2 - Ifges(5,1) + Ifges(5,2)) * t3) * t3 + ((2 * mrSges(6,3) + t6) * pkin(7)) + 0.2e1 * t5 * pkin(2); mrSges(3,1) + t5; mrSges(4,1) * t2 + mrSges(4,2) * t4 + mrSges(3,2); m(4); ((m(6) * pkin(4) + mrSges(5,1)) * t3 - mrSges(5,2) * t1) / t3; mrSges(5,3) + mrSges(6,3) + t6; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t7;
