% Return the minimum parameter vector for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t79 = -pkin(7) * m(6) - mrSges(6,3);
t71 = (m(5) + m(6));
t73 = (pkin(7) ^ 2);
t75 = (pkin(4) ^ 2);
t78 = (Ifges(5,2) + (t73 + t75) * m(6));
t77 = (mrSges(5,3) - t79);
t76 = 2 * pkin(7) * mrSges(6,3) + 2 * pkin(6) * t77 + Ifges(6,2) + t78;
t74 = pkin(6) ^ 2;
t70 = cos(pkin(8));
t69 = sin(pkin(8));
t1 = [Ifges(2,3) + Ifges(3,1) + t70 ^ 2 * (t74 * t71 + Ifges(4,1) + t76) + (-0.2e1 * t70 * Ifges(4,4) + (Ifges(4,2) + (pkin(3) ^ 2 + t74) * t71 + t76) * t69) * t69; mrSges(2,1); mrSges(2,2); mrSges(3,2); mrSges(3,3); m(3); pkin(3) * t71 + mrSges(4,1); mrSges(4,2); pkin(6) * t71 + mrSges(4,3) + t77; m(4) + t71; m(6) * t73 + Ifges(5,1) - t78; Ifges(5,4); t79 * pkin(4) + Ifges(5,5); Ifges(5,6); t75 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
