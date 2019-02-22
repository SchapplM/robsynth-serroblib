% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:29:14
% EndTime: 2019-02-22 11:29:14
% DurationCPUTime: 0.14s
% Computational Cost: add. (128->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t142 = cos(pkin(6));
t147 = cos(qJ(2));
t148 = cos(qJ(1));
t149 = t148 * t147;
t144 = sin(qJ(2));
t145 = sin(qJ(1));
t152 = t145 * t144;
t132 = -t142 * t149 + t152;
t140 = qJ(4) + pkin(11);
t138 = sin(t140);
t139 = cos(t140);
t141 = sin(pkin(6));
t155 = t141 * t148;
t127 = -t132 * t138 + t139 * t155;
t150 = t148 * t144;
t151 = t145 * t147;
t133 = t142 * t150 + t151;
t143 = sin(qJ(6));
t146 = cos(qJ(6));
t163 = t127 * t143 + t133 * t146;
t162 = t127 * t146 - t133 * t143;
t159 = t138 * t143;
t158 = t138 * t146;
t157 = t141 * t145;
t156 = t141 * t147;
t154 = t143 * t144;
t153 = t144 * t146;
t126 = t132 * t139 + t138 * t155;
t135 = -t142 * t152 + t149;
t134 = t142 * t151 + t150;
t131 = -t138 * t156 + t142 * t139;
t130 = -t142 * t138 - t139 * t156;
t125 = t134 * t138 + t139 * t157;
t124 = -t134 * t139 + t138 * t157;
t123 = t125 * t146 + t135 * t143;
t122 = -t125 * t143 + t135 * t146;
t1 = [t162, -t134 * t143 + t135 * t158, 0, -t124 * t146, 0, t122; t123, -t132 * t143 + t133 * t158, 0, t126 * t146, 0, t163; 0 (t138 * t153 + t143 * t147) * t141, 0, t130 * t146, 0, -t131 * t143 + t141 * t153; -t163, -t134 * t146 - t135 * t159, 0, t124 * t143, 0, -t123; t122, -t132 * t146 - t133 * t159, 0, -t126 * t143, 0, t162; 0 (-t138 * t154 + t146 * t147) * t141, 0, -t130 * t143, 0, -t131 * t146 - t141 * t154; t126, -t135 * t139, 0, t125, 0, 0; t124, -t133 * t139, 0, -t127, 0, 0; 0, -t141 * t144 * t139, 0, t131, 0, 0;];
JR_rot  = t1;
