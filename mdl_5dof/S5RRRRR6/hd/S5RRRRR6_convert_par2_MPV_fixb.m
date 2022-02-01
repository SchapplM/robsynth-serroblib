% Return the minimum parameter vector for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t93 = -pkin(9) * m(6) - mrSges(6,3);
t79 = (mrSges(5,3) - t93);
t84 = (m(5) + m(6));
t92 = -pkin(8) * t84 - t79;
t83 = (m(4) + t84);
t90 = (mrSges(4,3) - t92);
t89 = (pkin(3) ^ 2);
t88 = (pkin(4) ^ 2);
t87 = (pkin(8) ^ 2);
t86 = (pkin(9) ^ 2);
t82 = (t87 + t89);
t81 = (t86 + t88);
t80 = m(3) + t83;
t1 = [pkin(1) ^ 2 * t80 + Ifges(2,3); pkin(1) * t80 + mrSges(2,1); mrSges(2,2); Ifges(3,3) + Ifges(4,2) + Ifges(5,2) + Ifges(6,2) + 2 * pkin(9) * mrSges(6,3) + t81 * m(6) + 2 * pkin(8) * t79 + t82 * t84 + 2 * pkin(7) * t90 + (pkin(2) ^ 2 + pkin(7) ^ 2) * t83; pkin(2) * t83 + mrSges(3,1); -pkin(7) * t83 + mrSges(3,2) - t90; Ifges(4,1) - Ifges(4,2) + (-t82 + t87) * t84; Ifges(4,4); t92 * pkin(3) + Ifges(4,5); Ifges(4,6); t89 * t84 + Ifges(4,3); pkin(3) * t84 + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) + (-t81 + t86) * m(6); Ifges(5,4); t93 * pkin(4) + Ifges(5,5); Ifges(5,6); t88 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
