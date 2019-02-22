% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:27:29
% EndTime: 2019-02-22 11:27:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (163->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->38)
t146 = cos(pkin(6));
t147 = sin(qJ(2));
t150 = cos(qJ(1));
t152 = t150 * t147;
t148 = sin(qJ(1));
t149 = cos(qJ(2));
t153 = t148 * t149;
t134 = t146 * t152 + t153;
t144 = pkin(11) + qJ(4);
t140 = sin(t144);
t142 = cos(t144);
t145 = sin(pkin(6));
t155 = t145 * t150;
t128 = -t134 * t142 + t140 * t155;
t151 = t150 * t149;
t154 = t148 * t147;
t133 = -t146 * t151 + t154;
t143 = pkin(12) + qJ(6);
t139 = sin(t143);
t141 = cos(t143);
t165 = t128 * t139 + t133 * t141;
t164 = t128 * t141 - t133 * t139;
t161 = t139 * t142;
t160 = t141 * t142;
t159 = t142 * t149;
t158 = t145 * t147;
t157 = t145 * t148;
t156 = t145 * t149;
t126 = -t134 * t140 - t142 * t155;
t136 = -t146 * t154 + t151;
t135 = t146 * t153 + t152;
t132 = t146 * t140 + t142 * t158;
t131 = -t140 * t158 + t146 * t142;
t130 = t136 * t142 + t140 * t157;
t129 = t136 * t140 - t142 * t157;
t125 = t130 * t141 + t135 * t139;
t124 = -t130 * t139 + t135 * t141;
t1 = [t164, -t135 * t160 + t136 * t139, 0, -t129 * t141, 0, t124; t125, -t133 * t160 + t134 * t139, 0, t126 * t141, 0, t165; 0 (t139 * t147 + t141 * t159) * t145, 0, t131 * t141, 0, -t132 * t139 - t141 * t156; -t165, t135 * t161 + t136 * t141, 0, t129 * t139, 0, -t125; t124, t133 * t161 + t134 * t141, 0, -t126 * t139, 0, t164; 0 (-t139 * t159 + t141 * t147) * t145, 0, -t131 * t139, 0, -t132 * t141 + t139 * t156; t126, -t135 * t140, 0, t130, 0, 0; t129, -t133 * t140, 0, -t128, 0, 0; 0, t140 * t156, 0, t132, 0, 0;];
JR_rot  = t1;
