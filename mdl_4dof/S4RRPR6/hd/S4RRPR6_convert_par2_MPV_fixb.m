% Return the minimum parameter vector for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t72 = (pkin(6) ^ 2);
t73 = (pkin(3) ^ 2);
t78 = 2 * pkin(6) * mrSges(5,3) + Ifges(5,2);
t63 = Ifges(4,2) + (t72 + t73) * m(5) + t78;
t64 = m(5) * t72 + Ifges(4,1) + t78;
t80 = -t63 + t64;
t70 = sin(pkin(7));
t71 = cos(pkin(7));
t79 = t70 * t71;
t67 = t70 ^ 2;
t68 = t71 ^ 2;
t77 = t68 - t67;
t76 = Ifges(4,4) * t79;
t75 = -pkin(6) * m(5) - mrSges(5,3);
t66 = m(5) * pkin(3) + mrSges(4,1);
t74 = -t70 * mrSges(4,2) + t71 * t66;
t65 = pkin(3) * t75 + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t67 * t64 + 0.2e1 * t76 + t68 * t63 + (2 * pkin(5) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(5) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(5) * m(3) + mrSges(2,2) - mrSges(3,3); t80 * t77 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t76; t77 * Ifges(4,4) + t80 * t79 + Ifges(3,4); -t70 * Ifges(4,6) + t71 * t65 + Ifges(3,5); t71 * Ifges(4,6) + t70 * t65 + Ifges(3,6); (t73 * m(5)) + 0.2e1 * pkin(2) * t74 + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t74; t71 * mrSges(4,2) + t70 * t66 + mrSges(3,2); mrSges(4,3) - t75; m(4) + m(5); Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t1;
