% Return the minimum parameter vector for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t89 = -pkin(7) * m(6) - mrSges(6,3);
t80 = (m(5) + m(6));
t82 = (pkin(7) ^ 2);
t84 = (pkin(4) ^ 2);
t88 = (Ifges(5,2) + (t82 + t84) * m(6));
t87 = (mrSges(5,3) - t89);
t86 = 2 * pkin(7) * mrSges(6,3) + 2 * pkin(6) * t87 + Ifges(6,2) + t88;
t77 = sin(pkin(8));
t79 = cos(pkin(8));
t85 = t79 * mrSges(3,1) - t77 * mrSges(3,2);
t83 = pkin(6) ^ 2;
t78 = cos(pkin(9));
t76 = sin(pkin(9));
t1 = [Ifges(2,3) + Ifges(3,3) + t78 ^ 2 * (Ifges(4,2) + (pkin(3) ^ 2 + t83) * t80 + t86) + (0.2e1 * t78 * Ifges(4,4) + (t83 * t80 + Ifges(4,1) + t86) * t76) * t76 + 0.2e1 * t85 * pkin(1); mrSges(2,1) + t85; t77 * mrSges(3,1) + t79 * mrSges(3,2) + mrSges(2,2); m(3); pkin(3) * t80 + mrSges(4,1); mrSges(4,2); pkin(6) * t80 + mrSges(4,3) + t87; m(4) + t80; m(6) * t82 + Ifges(5,1) - t88; Ifges(5,4); t89 * pkin(4) + Ifges(5,5); Ifges(5,6); t84 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
