% Return the minimum parameter vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% MPV [19x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t80 = (pkin(6) ^ 2);
t86 = 2 * pkin(6) * mrSges(5,3) + Ifges(5,2);
t70 = m(5) * t80 + Ifges(4,1) + t86;
t81 = (pkin(3) ^ 2);
t73 = m(5) * t81 + Ifges(4,2);
t88 = t70 - t73;
t78 = sin(pkin(7));
t79 = cos(pkin(7));
t87 = t78 * t79;
t75 = t78 ^ 2;
t76 = t79 ^ 2;
t85 = t76 - t75;
t83 = pkin(6) * m(5) + mrSges(5,3);
t71 = t83 * pkin(3) + Ifges(4,4);
t84 = t71 * t87;
t72 = mrSges(4,2) - t83;
t74 = m(5) * pkin(3) + mrSges(4,1);
t82 = -t72 * t78 + t74 * t79;
t1 = [Ifges(2,3) + Ifges(3,2) + t75 * t70 + 0.2e1 * t84 + t76 * t73 + (2 * pkin(5) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(5) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(5) * m(3) + mrSges(2,2) - mrSges(3,3); t88 * t85 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t84; t85 * t71 + t88 * t87 + Ifges(3,4); Ifges(4,5) * t79 - Ifges(4,6) * t78 + Ifges(3,5); Ifges(4,5) * t78 + Ifges(4,6) * t79 + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t80 + t81) * m(5)) + 0.2e1 * t82 * pkin(2) + t86; mrSges(3,1) + t82; t72 * t79 + t74 * t78 + mrSges(3,2); mrSges(4,3); m(4) + m(5); Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t1;
