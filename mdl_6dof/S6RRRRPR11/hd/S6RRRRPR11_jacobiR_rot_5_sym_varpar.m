% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:23
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:23:23
% EndTime: 2019-02-22 12:23:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (112->31), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
t137 = cos(pkin(6));
t139 = sin(qJ(2));
t143 = cos(qJ(1));
t145 = t143 * t139;
t140 = sin(qJ(1));
t142 = cos(qJ(2));
t147 = t140 * t142;
t128 = t137 * t145 + t147;
t138 = sin(qJ(3));
t141 = cos(qJ(3));
t136 = sin(pkin(6));
t149 = t136 * t143;
t122 = -t128 * t141 + t138 * t149;
t144 = t143 * t142;
t148 = t140 * t139;
t127 = -t137 * t144 + t148;
t135 = qJ(4) + pkin(12);
t133 = sin(t135);
t134 = cos(t135);
t158 = t122 * t133 + t127 * t134;
t157 = t122 * t134 - t127 * t133;
t154 = t133 * t141;
t153 = t134 * t141;
t152 = t136 * t138;
t151 = t136 * t141;
t150 = t136 * t142;
t146 = t141 * t142;
t120 = -t128 * t138 - t141 * t149;
t130 = -t137 * t148 + t144;
t129 = t137 * t147 + t145;
t126 = t137 * t138 + t139 * t151;
t125 = t137 * t141 - t139 * t152;
t124 = t130 * t141 + t140 * t152;
t123 = t130 * t138 - t140 * t151;
t119 = t124 * t134 + t129 * t133;
t118 = -t124 * t133 + t129 * t134;
t1 = [t157, -t129 * t153 + t130 * t133, -t123 * t134, t118, 0, 0; t119, -t127 * t153 + t128 * t133, t120 * t134, t158, 0, 0; 0 (t133 * t139 + t134 * t146) * t136, t125 * t134, -t126 * t133 - t134 * t150, 0, 0; -t158, t129 * t154 + t130 * t134, t123 * t133, -t119, 0, 0; t118, t127 * t154 + t128 * t134, -t120 * t133, t157, 0, 0; 0 (-t133 * t146 + t134 * t139) * t136, -t125 * t133, -t126 * t134 + t133 * t150, 0, 0; t120, -t129 * t138, t124, 0, 0, 0; t123, -t127 * t138, -t122, 0, 0, 0; 0, t138 * t150, t126, 0, 0, 0;];
JR_rot  = t1;
