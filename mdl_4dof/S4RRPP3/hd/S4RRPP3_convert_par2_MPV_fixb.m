% Return the minimum parameter vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% MPV [16x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRPP3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t62 = sin(pkin(6));
t60 = t62 ^ 2;
t63 = cos(pkin(6));
t61 = t63 ^ 2;
t73 = t61 - t60;
t72 = t62 * t63;
t65 = Ifges(4,2) + Ifges(5,3);
t68 = Ifges(4,1) + Ifges(5,1);
t71 = t65 - t68;
t67 = Ifges(4,4) - Ifges(5,5);
t70 = t67 * t72;
t69 = t63 * mrSges(4,1) - t62 * mrSges(4,2);
t66 = Ifges(4,5) + Ifges(5,4);
t64 = Ifges(4,6) - Ifges(5,6);
t1 = [Ifges(2,3) + Ifges(3,2) + t60 * t68 + 0.2e1 * t70 + t61 * t65 + (2 * pkin(5) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(5) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(5) * m(3) + mrSges(2,2) - mrSges(3,3); -t73 * t71 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t70; t73 * t67 - t71 * t72 + Ifges(3,4); -t62 * t64 + t63 * t66 + Ifges(3,5); t62 * t66 + t63 * t64 + Ifges(3,6); 0.2e1 * pkin(2) * t69 + Ifges(5,2) + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t69; t62 * mrSges(4,1) + t63 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5);];
MPV = t1;
