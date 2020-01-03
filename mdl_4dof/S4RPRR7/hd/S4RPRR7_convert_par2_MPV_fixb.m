% Return the minimum parameter vector for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RPRR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t70 = (m(4) + m(5));
t74 = (pkin(3) ^ 2);
t78 = (m(5) * t74 + Ifges(4,2));
t77 = 2 * pkin(6) * mrSges(5,3) + Ifges(5,2);
t76 = 2 * pkin(5) * mrSges(4,3) + t78;
t75 = m(5) * pkin(6) + mrSges(5,3);
t73 = pkin(5) ^ 2;
t72 = pkin(6) ^ 2;
t69 = cos(pkin(7));
t68 = sin(pkin(7));
t1 = [Ifges(2,3) + t69 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t73) * t70 + t76) + (0.2e1 * t69 * Ifges(3,4) + (t70 * t73 + Ifges(3,1) + t76) * t68) * t68; mrSges(2,1); mrSges(2,2); pkin(2) * t70 + mrSges(3,1); mrSges(3,2); pkin(5) * t70 + mrSges(3,3) + mrSges(4,3); m(3) + t70; m(5) * t72 + Ifges(4,1) + t77 - t78; pkin(3) * t75 + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t72 + t74) * m(5) + t77; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t75; Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t1;
