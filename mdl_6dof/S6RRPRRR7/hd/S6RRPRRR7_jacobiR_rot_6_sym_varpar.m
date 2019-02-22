% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:43
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:43:25
% EndTime: 2019-02-22 11:43:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (98->25), mult. (164->24), div. (0->0), fcn. (248->8), ass. (0->28)
t138 = sin(qJ(4));
t140 = cos(qJ(2));
t155 = sin(qJ(2));
t156 = cos(qJ(4));
t129 = -t140 * t138 + t155 * t156;
t128 = t155 * t138 + t140 * t156;
t139 = sin(qJ(1));
t124 = t129 * t139;
t137 = qJ(5) + qJ(6);
t135 = sin(t137);
t154 = t124 * t135;
t136 = cos(t137);
t153 = t124 * t136;
t141 = cos(qJ(1));
t127 = t129 * t141;
t152 = t127 * t135;
t151 = t127 * t136;
t150 = t128 * t135;
t149 = t128 * t136;
t148 = t129 * t135;
t147 = t129 * t136;
t125 = t128 * t139;
t121 = -t125 * t136 - t141 * t135;
t142 = t125 * t135 - t141 * t136;
t126 = t128 * t141;
t123 = t126 * t136 - t139 * t135;
t122 = -t126 * t135 - t139 * t136;
t1 = [t121, -t151, 0, t151, t122, t122; t123, -t153, 0, t153, -t142, -t142; 0, t149, 0, -t149, -t148, -t148; t142, t152, 0, -t152, -t123, -t123; t122, t154, 0, -t154, t121, t121; 0, -t150, 0, t150, -t147, -t147; t124, -t126, 0, t126, 0, 0; -t127, -t125, 0, t125, 0, 0; 0, -t129, 0, t129, 0, 0;];
JR_rot  = t1;
