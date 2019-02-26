% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (128->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t141 = cos(pkin(6));
t146 = cos(qJ(2));
t147 = cos(qJ(1));
t148 = t147 * t146;
t143 = sin(qJ(2));
t144 = sin(qJ(1));
t151 = t144 * t143;
t131 = -t141 * t148 + t151;
t139 = pkin(11) + qJ(5);
t137 = sin(t139);
t138 = cos(t139);
t140 = sin(pkin(6));
t154 = t140 * t147;
t126 = -t131 * t137 + t138 * t154;
t149 = t147 * t143;
t150 = t144 * t146;
t132 = t141 * t149 + t150;
t142 = sin(qJ(6));
t145 = cos(qJ(6));
t162 = t126 * t142 + t132 * t145;
t161 = t126 * t145 - t132 * t142;
t158 = t137 * t142;
t157 = t137 * t145;
t156 = t140 * t144;
t155 = t140 * t146;
t153 = t142 * t143;
t152 = t143 * t145;
t125 = t131 * t138 + t137 * t154;
t134 = -t141 * t151 + t148;
t133 = t141 * t150 + t149;
t130 = -t137 * t155 + t141 * t138;
t129 = -t141 * t137 - t138 * t155;
t124 = t133 * t137 + t138 * t156;
t123 = -t133 * t138 + t137 * t156;
t122 = t124 * t145 + t134 * t142;
t121 = -t124 * t142 + t134 * t145;
t1 = [t161, -t133 * t142 + t134 * t157, 0, 0, -t123 * t145, t121; t122, -t131 * t142 + t132 * t157, 0, 0, t125 * t145, t162; 0 (t137 * t152 + t142 * t146) * t140, 0, 0, t129 * t145, -t130 * t142 + t140 * t152; -t162, -t133 * t145 - t134 * t158, 0, 0, t123 * t142, -t122; t121, -t131 * t145 - t132 * t158, 0, 0, -t125 * t142, t161; 0 (-t137 * t153 + t145 * t146) * t140, 0, 0, -t129 * t142, -t130 * t145 - t140 * t153; t125, -t134 * t138, 0, 0, t124, 0; t123, -t132 * t138, 0, 0, -t126, 0; 0, -t140 * t143 * t138, 0, 0, t130, 0;];
JR_rot  = t1;
