% Return the minimum parameter vector for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t73 = (m(4) + m(5));
t81 = 2 * pkin(7) * mrSges(5,3) + Ifges(5,2);
t80 = pkin(7) * m(5) + mrSges(5,3);
t79 = -pkin(6) * t73 - mrSges(4,3);
t78 = (mrSges(3,3) - t79);
t77 = (pkin(2) ^ 2);
t76 = (pkin(3) ^ 2);
t75 = (pkin(6) ^ 2);
t74 = pkin(7) ^ 2;
t71 = (m(3) + t73);
t70 = (t75 + t77);
t1 = [Ifges(2,3) + Ifges(3,2) + t76 * m(5) + Ifges(4,2) + 2 * pkin(6) * mrSges(4,3) + t70 * t73 + 2 * pkin(5) * t78 + (pkin(1) ^ 2 + pkin(5) ^ 2) * t71; pkin(1) * t71 + mrSges(2,1); -pkin(5) * t71 + mrSges(2,2) - t78; Ifges(3,1) - Ifges(3,2) + (-t70 + t75) * t73; Ifges(3,4); t79 * pkin(2) + Ifges(3,5); Ifges(3,6); t77 * t73 + Ifges(3,3); pkin(2) * t73 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (t74 - t76) * m(5) + t81; t80 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t74 + t76) * m(5) + t81; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t80; Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t1;
