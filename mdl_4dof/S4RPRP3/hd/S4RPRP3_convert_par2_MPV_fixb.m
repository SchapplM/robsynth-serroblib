% Return the minimum parameter vector for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m [5x1]
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
% MPV [15x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RPRP3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t6 = -Ifges(4,2) - Ifges(5,2);
t1 = m(4) * pkin(5) - mrSges(3,2) + mrSges(4,3);
t2 = m(4) * pkin(2) + mrSges(3,1);
t3 = sin(pkin(6));
t4 = cos(pkin(6));
t5 = t1 * t3 + t2 * t4;
t7 = [((pkin(2) ^ 2 + pkin(5) ^ 2) * m(4)) + (2 * mrSges(4,3) * pkin(5)) + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t5 * pkin(1) - t6; mrSges(2,1) + t5; -t1 * t4 + t2 * t3 + mrSges(2,2); m(3) + m(4); Ifges(4,1) + Ifges(5,1) + t6; Ifges(4,4) + Ifges(5,4); Ifges(4,5) + Ifges(5,5); Ifges(4,6) + Ifges(5,6); Ifges(4,3) + Ifges(5,3); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5);];
MPV = t7;
