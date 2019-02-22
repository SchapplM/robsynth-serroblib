% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP13_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:38:36
% EndTime: 2019-02-22 11:38:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (77->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t137 = cos(pkin(6));
t144 = cos(qJ(2));
t145 = cos(qJ(1));
t146 = t145 * t144;
t140 = sin(qJ(2));
t141 = sin(qJ(1));
t149 = t141 * t140;
t130 = -t137 * t146 + t149;
t139 = sin(qJ(4));
t143 = cos(qJ(4));
t136 = sin(pkin(6));
t154 = t136 * t145;
t125 = -t130 * t139 + t143 * t154;
t147 = t145 * t140;
t148 = t141 * t144;
t131 = t137 * t147 + t148;
t138 = sin(qJ(5));
t142 = cos(qJ(5));
t160 = t125 * t138 + t131 * t142;
t159 = t125 * t142 - t131 * t138;
t156 = t136 * t143;
t155 = t136 * t144;
t153 = t138 * t139;
t152 = t138 * t140;
t151 = t139 * t142;
t150 = t140 * t142;
t124 = t130 * t143 + t139 * t154;
t133 = -t137 * t149 + t146;
t132 = t137 * t148 + t147;
t129 = t137 * t143 - t139 * t155;
t128 = -t137 * t139 - t143 * t155;
t123 = t132 * t139 + t141 * t156;
t122 = t141 * t136 * t139 - t132 * t143;
t121 = t123 * t142 + t133 * t138;
t120 = -t123 * t138 + t133 * t142;
t1 = [t159, -t132 * t138 + t133 * t151, 0, -t122 * t142, t120, 0; t121, -t130 * t138 + t131 * t151, 0, t124 * t142, t160, 0; 0 (t138 * t144 + t139 * t150) * t136, 0, t128 * t142, -t129 * t138 + t136 * t150, 0; -t160, -t132 * t142 - t133 * t153, 0, t122 * t138, -t121, 0; t120, -t130 * t142 - t131 * t153, 0, -t124 * t138, t159, 0; 0 (-t139 * t152 + t142 * t144) * t136, 0, -t128 * t138, -t129 * t142 - t136 * t152, 0; t124, -t133 * t143, 0, t123, 0, 0; t122, -t131 * t143, 0, -t125, 0, 0; 0, -t140 * t156, 0, t129, 0, 0;];
JR_rot  = t1;
